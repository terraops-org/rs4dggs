#!/usr/bin/env bash
# The local gate: a milestone is done only when this passes. Actions minutes cost money,
# so this, not CI, is where verification happens.
set -euo pipefail
cd "$(dirname "$0")/.."

# The crate reads two environment variables, each of which lowers one of its limits for the whole
# process, and the tests take both limits at their ceilings: the oracle suites assert it, and
# under a lowered limit on a list of sub-zones the unit tests that build lists fail with the
# crate's own refusal, which does not name the variable. The gate therefore refuses to run while
# either is set, whatever its value, rather than fail for a reason it would not state.
for var in RS4DGGS_MAX_MATERIALISED_SUB_ZONES RS4DGGS_MAX_EDGE_REFINEMENT; do
  [[ -z "${!var+set}" ]] || {
    echo "$var is set in the environment (to \"${!var}\")." >&2
    echo "It lowers a limit of the crate, and the tests assume that limit at its ceiling, so that" >&2
    echo "under it they may fail without naming the cause. Unset it and run the gate again:" >&2
    echo "  unset $var" >&2
    exit 1
  }
done

# `--without-dggal` checks everything save what needs DGGAL's own libraries: the build identity and
# the oracle suites. GitHub's workflows run this form, since the runners cannot install the DGGAL
# build the compiled operation orders were read from; every commit runs the full gate locally, in
# the pre-commit hook.
WITH_DGGAL=1
if [[ "${1:-}" == "--without-dggal" ]]; then WITH_DGGAL=0; fi
if (( WITH_DGGAL )); then
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
else
  echo "== DGGAL: not used (--without-dggal); the build identity and the oracle suites are skipped"
fi

# Our crates only: the vendored DGGAL bindings under crates/dggal-oracle/vendor are unmodified
# upstream code and must not be reformatted, linted or doctested. -p already scopes fmt and
# test to our own crates, but Cargo runs clippy's RUSTC_WORKSPACE_WRAPPER over every
# workspace member regardless of -p, so without --no-deps the vendored path dependencies
# would be linted (and fail) too; --no-deps keeps clippy's lint pass on the selected
# packages only, while the vendored crates still build normally for linking.
CRATES=(-p rs4dggs -p rs4dggs-cli -p rs4dggs-ogc -p dggal-oracle)
FEATURES=(--features rs4dggs/oracle,rs4dggs-cli/oracle,rs4dggs-ogc/oracle)
if (( ! WITH_DGGAL )); then CRATES=(-p rs4dggs -p rs4dggs-cli -p rs4dggs-ogc); FEATURES=(); fi
# `cargo fmt` is scoped to the owned crates deliberately. The four crates under
# crates/dggal-oracle/vendor are verbatim upstream DGGAL and eCere binding sources, kept
# byte-identical so that they can be diffed against a new release; reformatting them would
# destroy that property and create a diff nobody can review.
echo "== fmt";     cargo fmt "${CRATES[@]}" --check
# --features rs4dggs/oracle, not a bare --features oracle: with several packages selected by
# -p, an unqualified feature name is only accepted when every selected package declares
# it, and the others do not. Each crate's oracle feature is its own; the qualified form
# says so. In rs4dggs it reaches the inline oracle-dependent tests, in grid.rs,
# topologies/hex_a7.rs, topologies/hex_a3.rs and topologies/hex_a3_subzones.rs, and the six
# integration suites under tests/, one per grid; in rs4dggs-cli and rs4dggs-ogc, the suite that
# compares each with DGGAL's own `dgg`,
# so that clippy and test cover them here, in the one place the feature is meaningful,
# exactly as they did before it existed.
echo "== clippy";  cargo clippy "${CRATES[@]}" --all-targets --no-deps "${FEATURES[@]}" -- -D warnings
echo "== test";    cargo test "${CRATES[@]}" "${FEATURES[@]}"
echo "== rs4dggs-ogc: definitions and goldens"
# Every JSON and GeoJSON file the encodings keep, the definitions of the grids and the documents
# the tests compare with, must be JSON as RFC 8259 has it, in UTF-8: Python's own reader accepts
# `NaN` and `Infinity`, which are not JSON, and keeps the last of two members of one name in
# silence, so that both are refused here. A directory not yet present holds no file.
python3 - crates/rs4dggs-ogc/definitions crates/rs4dggs-ogc/tests/goldens <<'EOF'
import json, pathlib, sys

def no_constant(token):
    raise ValueError(f"{token} is not JSON")

def no_repeated_member(pairs):
    names = [k for k, _ in pairs]
    if len(set(names)) != len(names):
        raise ValueError(f"a member is repeated among {names}")
    return dict(pairs)

for d in map(pathlib.Path, sys.argv[1:]):
    files = sorted(p for p in d.rglob("*") if p.suffix in (".json", ".geojson")) if d.is_dir() else []
    for p in files:
        try:
            json.loads(p.read_bytes().decode("utf-8"), parse_constant=no_constant,
                       object_pairs_hook=no_repeated_member)
        except ValueError as e:
            sys.exit(f"{p}: {e}")
    print(f"{d}: {len(files)} files parsed")
EOF
# The goldens of the zone lists, at least: a directory emptied by mistake would pass the step above.
goldens=$(find crates/rs4dggs-ogc/tests/goldens -type f -name '*.json' | wc -l)
if (( goldens < 2 )); then
  echo "crates/rs4dggs-ogc/tests/goldens: $goldens JSON goldens, fewer than the 2 of the zone lists" >&2; exit 1
fi
# The six definitions of the grids, and the two that are the OGC register's files, held to the
# register's bytes, so that a verbatim copy cannot be edited unawares.
definitions=$(find crates/rs4dggs-ogc/definitions -type f -name '*.json' | wc -l)
if (( definitions < 6 )); then
  echo "crates/rs4dggs-ogc/definitions: $definitions definitions, fewer than the six grids" >&2; exit 1
fi
sha256sum --check --quiet <<'EOF'
bfc76a11c3a8b8f7c642f8efd997b664f4c03588ef91afbca8728bb5e9ec4bb2  crates/rs4dggs-ogc/definitions/ISEA3H.json
bba46a41a7fa19786ddd80081bbaea5b6da180d1e3b7de51686a4534eddfbc65  crates/rs4dggs-ogc/definitions/IVEA3H.json
ef345a7920dd0ac58a342ee8867e455c2930c8d0baa6fc64b924be2929538a42  crates/rs4dggs-ogc/tests/fixtures/ogc/dggs-json.json
d41a00f8fa473e231eaa93bd43dd72ecdf8bbfa69c66a9e77f2e974bca2a61af  crates/rs4dggs-ogc/tests/fixtures/ogc/1-temperature.json
EOF
echo "crates/rs4dggs-ogc/definitions: the two registered definitions are the register's bytes"
# The goldens of zone data, at least the eight of DGGS-JSON.
data_goldens=$(find crates/rs4dggs-ogc/tests/goldens -type f -name 'data-*.json' | wc -l)
if (( data_goldens < 8 )); then
  echo "crates/rs4dggs-ogc/tests/goldens: $data_goldens goldens of zone data, fewer than 8" >&2; exit 1
fi
# And the two of DGGS-UBJSON, which no step above reads.
ubjson_goldens=$(find crates/rs4dggs-ogc/tests/goldens -type f -name 'data-*.ubj' | wc -l)
if (( ubjson_goldens < 2 )); then
  echo "crates/rs4dggs-ogc/tests/goldens: $ubjson_goldens goldens of DGGS-UBJSON, fewer than 2" >&2; exit 1
fi
if (( WITH_DGGAL )); then
echo "== rs4dggs-ogc: the DGGS-JSON schema"
# Every golden of zone data, and the standard's own example, against the standard's schema of
# DGGS-JSON. Run with DGGAL alone, as the oracle suites are: it needs Python's `jsonschema`, which
# GitHub's runners do not have, and an absent module fails the gate rather than skip the step.
python3 - crates/rs4dggs-ogc/tests/fixtures/ogc crates/rs4dggs-ogc/tests/goldens <<'EOF'
import json, pathlib, sys

try:
    from jsonschema import Draft202012Validator
except ImportError:
    sys.exit("the DGGS-JSON schema step needs the Python module jsonschema: pip install jsonschema")

fixtures, goldens = map(pathlib.Path, sys.argv[1:])
def read(p):
    return json.loads(p.read_text(encoding="utf-8"))

published = Draft202012Validator(read(fixtures / "dggs-json.json"))
schema = read(fixtures / "dggs-json.json")
# The published schema marks the values `nullable`, a word of OpenAPI 3.0 that a validator of
# JSON Schema Draft 2020-12 ignores, so that it refuses the `null` the standard requires for a
# missing value. Patched here, in memory, in that one place; the file stays as published.
schema["properties"]["values"]["additionalProperties"]["items"]["properties"]["data"]["items"] = {
    "type": ["number", "null"]}
patched = Draft202012Validator(schema)

documents = sorted(goldens.glob("data-*.json"))
if not documents:
    sys.exit(f"{goldens}: no golden of zone data to validate")
for p in [fixtures / "1-temperature.json", *documents]:
    errors = [e.message for e in patched.iter_errors(read(p))]
    if errors:
        sys.exit(f"{p}: {errors}")

# The step can fail: a document without its zone is refused, and so is a `null` by the schema as
# published.
first = read(documents[0])
del first["zoneId"]
if patched.is_valid(first):
    sys.exit("the schema accepts a document without zoneId")
if published.is_valid(read(goldens / "data-ISEA3H-C2-23-C-depth1.null-at-2.json")):
    sys.exit("the published schema accepts a null value: the patch may no longer be needed")
print(f"{goldens}: the example and {len(documents)} goldens of zone data are valid DGGS-JSON")
EOF
fi
# The GeoJSON goldens of a zone, of its feature, of the zone lists and of zone data, at least.
geo_goldens=$(find crates/rs4dggs-ogc/tests/goldens -type f -name '*.geojson' | wc -l)
if (( geo_goldens < 7 )); then
  echo "crates/rs4dggs-ogc/tests/goldens: $geo_goldens GeoJSON goldens, fewer than the 7 of a zone, its lists and zone data" >&2; exit 1
fi
echo "== doc";      RUSTDOCFLAGS="-D warnings" cargo doc -p rs4dggs --no-deps
RUSTDOCFLAGS="-D warnings" cargo doc -p rs4dggs-ogc --no-deps
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
# The same holds for the encodings: the rs4dggs on crates.io lacks functions they call.
cargo package -p rs4dggs-ogc --allow-dirty --list >/dev/null
echo "rs4dggs-ogc: file set checked; its full package check waits for rs4dggs on crates.io"
echo "== README"
# The repository's front page and the crate's are one document kept in two places: GitHub
# renders the root copy, crates.io and docs.rs the crate's. The crate's is the one to edit.
if ! cmp -s README.md crates/rs4dggs/README.md; then
  echo "README.md differs from crates/rs4dggs/README.md; edit the crate's copy and copy it to the root" >&2
  exit 1
fi
echo "== without DGGAL"
# A clone without the engine must still build, test, document and run the library, the tool and
# the encodings: only the test-only oracle needs DGGAL. A target directory of its own, so that the
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
out=$(ex hierarchy -- 00); grep -q '11 children' <<<"$out"
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
echo "== declared minimum (rust-version = 1.85, rs4dggs, rs4dggs-cli and rs4dggs-ogc; see Cargo.toml)"
# Not the whole gate above: the test-only oracle (dggal-oracle, its vendored ecrt) needs
# rust 1.87 for integer_sign_cast, so this checks the published crates on their own, and
# only when 1.85 is actually installed. Never silent: absent, it says so rather than
# skipping without a trace.
if ! command -v rustup >/dev/null 2>&1; then
  echo "declared minimum not checked: rustup is not on the PATH, so no toolchain can be selected by name"
elif ! rustup toolchain list 2>/dev/null | grep -q '^1\.85'; then
  echo "declared minimum not checked: rust 1.85 is not installed (rustup toolchain install 1.85 --profile minimal)"
else
  cargo +1.85 check -p rs4dggs -p rs4dggs-cli -p rs4dggs-ogc
fi
# Local additions this checkout keeps outside the published set, run here, where present.
if [[ -x check.local.sh ]]; then ./check.local.sh; fi
echo "check.sh: all green"
