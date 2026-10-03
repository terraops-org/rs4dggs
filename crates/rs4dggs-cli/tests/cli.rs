//! The tool run in process, as a person would run it.

use std::io::Cursor;

/// Runs the tool with `args` and `stdin`, and answers its exit code, standard output and standard
/// error.
pub fn run(args: &[&str], stdin: &str) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
    let mut input = Cursor::new(stdin.as_bytes().to_vec());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = rs4dggs_cli::run(&args, &mut input, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

#[test]
fn no_arguments_prints_the_help() {
    let (code, out, _) = run(&[], "");
    assert_eq!(code, 0);
    assert!(out.contains("rs4dggs <grid> <command>"), "{out}");
}

#[test]
fn help_and_version() {
    for args in [&["help"][..], &["--help"], &["-h"], &["igeo7", "help"]] {
        let (code, out, _) = run(args, "");
        assert_eq!(code, 0, "{args:?}");
        assert!(out.contains("Commands"), "{args:?}: {out}");
    }
    let (code, out, _) = run(&["igeo7", "help", "zone"], "");
    assert_eq!(code, 0);
    assert!(out.contains("rs4dggs igeo7 zone 38.72,-9.14 10"), "{out}");
    let (code, out, _) = run(&["--version"], "");
    assert_eq!(code, 0);
    let version = format!("rs4dggs-cli {}\n", env!("CARGO_PKG_VERSION"));
    assert_eq!(out, version);
}

#[test]
fn grids_lists_the_six() {
    let (code, out, _) = run(&["grids"], "");
    assert_eq!(code, 0);
    for name in ["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"] {
        assert!(out.contains(name), "{name}: {out}");
    }
    assert!(out.contains("0 to 19") && out.contains("0 to 33"), "{out}");
}

#[test]
fn grid_info_and_aliases() {
    for spelling in ["igeo7", "IGEO7", "isea7h_z7", "ISEA7H_Z7"] {
        let (code, out, _) = run(&[spelling, "info"], "");
        assert_eq!(code, 0, "{spelling}");
        assert!(out.starts_with("IGEO7\n"), "{spelling}: {out}");
        assert!(out.contains("aperture     7"), "{out}");
    }
    let (code, out, _) = run(&["ivea7h_z7"], "");
    assert_eq!(code, 0);
    assert!(out.starts_with("IVEA7H\n"), "{out}");
}

#[test]
fn usage_errors_exit_2_and_say_what_to_do() {
    let (code, _, err) = run(&["h3", "info"], "");
    assert_eq!(code, 2);
    assert!(
        err.contains("unknown grid") && err.contains("igeo7"),
        "{err}"
    );
    let (code, _, err) = run(&["igeo7", "frobnicate"], "");
    assert_eq!(code, 2);
    assert!(
        err.contains("unknown command") && err.contains("rs4dggs help"),
        "{err}"
    );
    let (code, _, err) = run(&["igeo7", "zone"], "");
    assert_eq!(code, 2);
    assert!(err.contains("38.72,-9.14"), "{err}");
    let (code, _, err) = run(&["igeo7", "info", "-f", "yaml"], "");
    assert_eq!(code, 2);
    assert!(err.contains("text, csv or geojson"), "{err}");
    let (code, _, err) = run(&["igeo7", "info", "--frobnicate"], "");
    assert_eq!(code, 2);
    assert!(err.contains("unknown option"), "{err}");
}

#[test]
fn a_negative_number_is_an_argument_not_an_option() {
    let inv = rs4dggs_cli::parse_for_tests(&["igeo7", "zone", "38.72", "-9.14", "10"]).unwrap();
    assert_eq!(inv, "zone 38.72,-9.14 res Some(10)");
    let inv = rs4dggs_cli::parse_for_tests(&["igeo7", "zone", "-38.72,-9.14"]).unwrap();
    assert_eq!(inv, "zone -38.72,-9.14 res None");
    let inv = rs4dggs_cli::parse_for_tests(&["igeo7", "zone", "-", "--format", "csv"]).unwrap();
    assert_eq!(inv, "zone - res None");
}

#[test]
fn zone_card_on_igeo7() {
    let (code, out, err) = run(&["igeo7", "zone", "38.72,-9.14", "10"], "");
    assert_eq!(code, 0, "{err}");
    assert!(
        out.starts_with(
            "IGEO7 zone 006415654636  (integer 940643638281502719, hex 0D0DD667BFFFFFFF)\n"
        ),
        "{out}"
    );
    assert!(out.contains("  resolution  10, hexagon\n"), "{out}");
    assert!(
        out.contains("  parents     2: 00641565463, 00641565462\n"),
        "{out}"
    );
    assert!(
        out.contains("  children    13: 0064156546360 (centre), 0064156546361, 0064156546365, 0064156546364, 0064156546366, 0064156546362, 0064156546363, 0064156546342, 0064156546033, 0064156546251, 0064156546215, 0064156546324, 0064156546306\n"),
        "{out}"
    );
    assert!(out.contains("  neighbours  6: 006415654634, 006415654603, 006415654625, 006415654621, 006415654632, 006415654630\n"), "{out}");
    assert!(out.contains("  vertices    6:\n"), "{out}");
    // The same card from the zone's text, its decimal and its hexadecimal identifiers.
    for z in [
        "006415654636",
        "940643638281502719",
        "0x0D0DD667BFFFFFFF",
        "0xd0dd667bfffffff",
    ] {
        let (code, again, err) = run(&["igeo7", "info", z], "");
        assert_eq!((code, again.as_str()), (0, out.as_str()), "{z}: {err}");
    }
}

#[test]
fn zone_card_on_isea3h_marks_the_centroid_parent() {
    let (code, out, err) = run(&["isea3h", "zone", "38.72", "-9.14", "5"], "");
    assert_eq!(code, 0, "{err}");
    assert!(
        out.starts_with(
            "ISEA3H zone C2-23-C  (integer 306244774661193870, hex 044000000000008E)\n"
        ),
        "{out}"
    );
    assert!(
        out.contains("  parents     3: C2-23-A, C4-5-A, C4-6-A (centroid parent)\n"),
        "{out}"
    );
    assert!(out.contains("  children    7: D2-128-A (centre), D2-10C-A, D2-10D-A, D4-11-A, D4-10-A, D2-143-A, D2-127-A\n"), "{out}");
}

#[test]
fn every_resolution_without_one() {
    let (code, out, _) = run(&["igeo7", "zone", "38.72,-9.14"], "");
    assert_eq!(code, 0);
    assert!(out.contains("\n    5  0064156\n"), "{out}");
    assert!(out.contains("\n   10  006415654636\n"), "{out}");
    assert_eq!(
        out.lines()
            .filter(|l| l.trim_start().starts_with(char::is_numeric))
            .count(),
        20,
        "{out}"
    );
}

#[test]
fn precision_rounds_every_coordinate() {
    let (_, out, _) = run(&["igeo7", "info", "006415654636", "-precision", "6"], "");
    assert!(
        out.contains("  centroid    38.719638, -9.140392   (lat, lon)\n"),
        "{out}"
    );
}

#[test]
fn input_errors_exit_1() {
    for (args, says) in [
        (&["igeo7", "zone", "95,0", "3"][..], "latitude"),
        (&["igeo7", "zone", "0,190", "3"], "longitude"),
        (&["igeo7", "zone", "0,0", "20"], "19"),
        (&["igeo7", "info", "(null)"], "null zone"),
        (&["igeo7", "info", "C2-23-C"], "IGEO7"),
        (&["isea3h", "info", "Z9-0-A"], "ISEA3H"),
        (&["igeo7", "info", "0064156546360000000000"], "IGEO7"),
        (&["igeo7", "info", "0xZZ"], "hexadecimal"),
        (&["isea3h", "sub", "A4-0-A", "-depth", "20"], "-depth"),
    ] {
        let (code, _, err) = run(args, "");
        assert_eq!(code, 1, "{args:?}: {err}");
        assert!(err.contains(says), "{args:?}: {err}");
    }
}

#[test]
fn neighbours_disk_rel_sub_index() {
    let (_, out, _) = run(&["igeo7", "neighbours", "006415654636"], "");
    assert_eq!(
        out,
        "IGEO7 zone 006415654636: 6 neighbours\n  006415654634\n  006415654603\n  006415654625\n  006415654621\n  006415654632\n  006415654630\n"
    );
    let (_, out, _) = run(&["igeo7", "disk", "006415654636", "2"], "");
    assert!(out.starts_with("IGEO7 disk of radius 2 about 006415654636\n  ring 0 (1): 006415654636\n  ring 1 (6): 006415654634, "), "{out}");
    assert!(out.contains("\n  ring 2 (12): "), "{out}");
    let (_, out, _) = run(&["igeo7", "rel", "006415654636", "006415654634"], "");
    assert!(out.contains("are neighbours"), "{out}");
    let (_, out, _) = run(&["isea3h", "rel", "A4-0-A", "B2-5-C"], "");
    assert!(
        out.contains("A4-0-A is coarser than B2-5-C by 3 resolutions"),
        "{out}"
    );
    assert!(out.contains("A4-0-A is an ancestor of B2-5-C"), "{out}");
    assert!(
        out.contains("B2-5-C is sub-zone 8 of A4-0-A, at depth 3"),
        "{out}"
    );
    let (_, out, _) = run(&["isea3h", "sub", "A4-0-A", "-depth", "3"], "");
    assert!(
        out.starts_with(
            "ISEA3H zone A4-0-A: 31 sub-zones at depth 3\n      0  B2-7-B\n      1  B2-4-D\n"
        ),
        "{out}"
    );
    let (_, out, _) = run(&["isea3h", "sub", "A4-0-A", "8", "-depth", "3"], "");
    assert_eq!(out, "sub-zone 8 of A4-0-A at depth 3: B2-5-C\n");
    let (_, out, _) = run(&["isea3h", "index", "A4-0-A", "B2-5-C"], "");
    assert_eq!(out, "B2-5-C is sub-zone 8 of A4-0-A, at depth 3\n");
    let (code, out, _) = run(&["isea3h", "index", "A4-0-A", "B6-5-C"], "");
    assert_eq!(
        (code, out.as_str()),
        (0, "B6-5-C has no index among the sub-zones of A4-0-A\n")
    );
    // The aperture-7 grids order their sub-zones as well, in the engine's own scanlines.
    let (code, out, _) = run(&["igeo7", "sub", "006415654636"], "");
    assert_eq!(code, 0, "{out}");
    assert!(
        out.starts_with(
            "IGEO7 zone 006415654636: 13 sub-zones at depth 1\n      0  0064156546306\n      1  0064156546342\n"
        ),
        "{out}"
    );
    let (code, out, _) = run(&["igeo7", "index", "0064156", "00641565"], "");
    assert_eq!(
        (code, out.as_str()),
        (0, "00641565 is sub-zone 5 of 0064156, at depth 1\n")
    );
}

#[test]
fn csv_for_one_zone_and_in_batch() {
    let (code, out, _) = run(&["igeo7", "zone", "38.72,-9.14", "10", "-f", "csv"], "");
    assert_eq!(
        (code, out.as_str()),
        (0, "lat,lon,zone\n38.72,-9.14,006415654636\n")
    );
    let input = "lat,lon\n# Lisbon and Porto\n38.72,-9.14\r\n\n41.15 -8.61\nnonsense\n95,0\n";
    let (code, out, err) = run(&["igeo7", "zone", "-", "5"], input);
    assert_eq!(code, 1, "{err}");
    assert_eq!(
        out,
        "lat,lon,zone\n38.72,-9.14,0064156\n41.15,-8.61,0065442\n"
    );
    assert!(err.contains("line 6") && err.contains("line 7"), "{err}");
    assert!(!err.contains("line 1"), "the header is not an error: {err}");
    let (_, out, _) = run(&["igeo7", "neighbours", "-"], "006415654636");
    assert_eq!(out.lines().next(), Some("zone,neighbour"));
    assert_eq!(out.lines().count(), 7);
    let (_, out, _) = run(&["igeo7", "disk", "-", "1"], "006415654636\n");
    assert_eq!(out.lines().next(), Some("zone,ring,member"));
    assert_eq!(out.lines().count(), 8);
    let (_, out, _) = run(&["isea3h", "info", "-"], "C2-23-C\n");
    assert_eq!(
        out.lines().nth(1).unwrap().split(',').next_back(),
        Some("C2-23-A C4-5-A C4-6-A")
    );
}

#[test]
fn batch_reads_bytes_that_are_not_utf8() {
    let args: Vec<String> = ["igeo7", "zone", "-", "3"]
        .iter()
        .map(|a| a.to_string())
        .collect();
    let mut input = std::io::Cursor::new(b"38.72,-9.14\n\xff\xfe,1\n".to_vec());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(rs4dggs_cli::run(&args, &mut input, &mut out, &mut err), 1);
    assert_eq!(String::from_utf8(out).unwrap().lines().count(), 2);
}

#[test]
fn geojson_features() {
    let (code, out, _) = run(&["igeo7", "geom", "006415654636"], "");
    assert_eq!(code, 0);
    assert!(out.starts_with(r#"{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"zone":"006415654636","grid":"IGEO7","resolution":10},"geometry":{"type":"Polygon","coordinates":[[["#), "{out}");
    let (_, out, _) = run(
        &[
            "igeo7",
            "disk",
            "006415654636",
            "1",
            "-f",
            "geojson",
            "-centroids",
        ],
        "",
    );
    assert_eq!(out.matches(r#""type":"Point""#).count(), 7);
    assert!(out.contains(r#""ring":1"#));
    // Base cell 0's zone about the north pole, at resolution 5: one polygon reaching the pole.
    let (_, out, _) = run(&["igeo7", "zone", "90,0", "5", "-f", "geojson"], "");
    assert!(
        out.contains("[180,90]") && out.contains("[-180,90]"),
        "{out}"
    );
    let (code, _, err) = run(&["igeo7", "geom", "-", "-f", "csv"], "00\n");
    assert_eq!(code, 2, "{err}");
}

#[test]
fn output_to_a_file_and_to_nowhere() {
    let dir = std::env::temp_dir().join(format!("rs4dggs-cli-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("cell.geojson");
    let (code, out, _) = run(
        &[
            "igeo7",
            "geom",
            "006415654636",
            "-o",
            path.to_str().unwrap(),
        ],
        "",
    );
    assert_eq!((code, out.as_str()), (0, ""));
    assert!(
        std::fs::read_to_string(&path)
            .unwrap()
            .starts_with(r#"{"type":"FeatureCollection""#)
    );
    let (code, _, err) = run(&["igeo7", "info", "-o", "/nonexistent-dir/x.txt"], "");
    assert_eq!(code, 1);
    assert!(err.contains("cannot write"), "{err}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_closed_pipe_ends_quietly() {
    struct Closed;
    impl std::io::Write for Closed {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::ErrorKind::BrokenPipe.into())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let args: Vec<String> = ["igeo7", "disk", "00", "30"]
        .iter()
        .map(|a| a.to_string())
        .collect();
    let mut err = Vec::new();
    assert_eq!(
        rs4dggs_cli::run(
            &args,
            &mut std::io::Cursor::new(Vec::new()),
            &mut Closed,
            &mut err
        ),
        0
    );
    assert!(err.is_empty());
}

#[test]
fn a_failing_command_leaves_an_existing_output_file_alone() {
    let dir = std::env::temp_dir().join(format!("rs4dggs-cli-keep-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("keep.txt");
    let p = path.to_str().unwrap();
    for args in [
        &["igeo7", "help", "nope", "-o", p][..],
        &["igeo7", "info", "nonsense", "-o", p],
        &["igeo7", "zone", "95,0", "5", "-f", "csv", "-o", p],
        &["igeo7", "geom", "nonsense", "-o", p],
        &["igeo7", "geom", "-", "-f", "csv", "-o", p],
    ] {
        std::fs::write(&path, "precious\n").unwrap();
        let (code, _, _) = run(args, "");
        assert_ne!(code, 0, "{args:?}");
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "precious\n",
            "{args:?}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_count_of_one_is_singular() {
    let (_, out, _) = run(&["igeo7", "rel", "00", "006"], "");
    assert!(
        out.contains("00 is coarser than 006 by 1 resolution\n"),
        "{out}"
    );
    let (_, out, _) = run(&["igeo7", "rel", "006", "00"], "");
    assert!(
        out.contains("006 is finer than 00 by 1 resolution\n"),
        "{out}"
    );
}

#[test]
fn geojson_of_batch_input_is_always_a_closed_collection() {
    let (code, out, err) = run(&["igeo7", "info", "-", "-f", "geojson"], "006415654636\n");
    assert_eq!(code, 0, "{err}");
    assert!(out.ends_with("]}\n") && err.is_empty(), "{out}{err}");
    let (code, out, err) = run(
        &["igeo7", "zone", "-", "5", "-f", "geojson"],
        "nonsense\n1,1\n",
    );
    assert_eq!((code, out.matches("Feature\"").count()), (0, 1), "{err}");
    let (code, out, _) = run(&["igeo7", "zone", "-", "5", "-f", "geojson"], "95,0\n");
    assert_eq!(
        (code, out.as_str()),
        (1, "{\"type\":\"FeatureCollection\",\"features\":[]}\n")
    );
}

#[test]
fn text_in_batch_is_separated_by_blank_lines_and_other_formats_are_refused() {
    let (_, out, _) = run(
        &["igeo7", "neighbours", "-", "-f", "text"],
        "006415654636\n006415654634\n",
    );
    assert_eq!(out.matches("\n\nIGEO7 zone").count(), 1, "{out}");
    let (code, _, err) = run(&["igeo7", "rel", "00", "006", "-f", "csv"], "");
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("rel writes text"), "{err}");
    let (code, _, err) = run(&["igeo7", "geom", "00", "-f", "text"], "");
    assert_eq!(
        (code, err.contains("geom writes geojson")),
        (2, true),
        "{err}"
    );
}

#[test]
fn text_batch_writes_a_blank_line_only_between_answers() {
    let input = "95,0\n006415654636\nzzz\n006415654634\n";
    let (code, out, err) = run(&["igeo7", "info", "-", "-f", "text"], input);
    assert_eq!(code, 1, "{err}");
    assert!(!out.starts_with('\n'), "{out}");
    assert_eq!(out.matches("\n\n").count(), 1, "{out}");
    assert_eq!(out.matches("IGEO7 zone").count(), 2, "{out}");
    assert!(!out.ends_with("\n\n"), "{out}");
}

#[test]
fn the_null_zone_is_left_out_of_geojson_and_answered_in_csv() {
    let (code, out, err) = run(
        &["igeo7", "zone", "89.999999,-120", "19", "-f", "geojson"],
        "",
    );
    assert_eq!(code, 0, "{err}");
    assert_eq!(out, "{\"type\":\"FeatureCollection\",\"features\":[]}\n");
    assert!(
        err.contains("1 input lies in no cell and was left out"),
        "{err}"
    );
    let (_, out, _) = run(&["igeo7", "zone", "89.999999,-120", "19", "-f", "csv"], "");
    assert_eq!(out, "lat,lon,zone\n89.999999,-120,\n");
    let (code, out, err) = run(
        &["igeo7", "zone", "-", "19", "-f", "geojson"],
        "38.72,-9.14\n89.999999,-120\n41.15,-8.61\n",
    );
    assert_eq!((code, out.matches("Feature\"").count()), (0, 2), "{err}");
}

#[test]
fn commands_with_two_arguments_do_not_read_standard_input() {
    for args in [
        &["igeo7", "rel", "-", "00"][..],
        &["igeo7", "sub", "-"],
        &["igeo7", "index", "-", "00"],
    ] {
        let (code, _, err) = run(args, "");
        assert_eq!(code, 2, "{args:?}");
        assert!(
            err.contains("does not read standard input"),
            "{args:?}: {err}"
        );
    }
}

#[test]
fn a_zone_with_every_vertex_on_a_pole_has_no_area_to_draw() {
    let (code, out, err) = run(&["igeo7", "zone", "90,0", "19", "-f", "geojson"], "");
    assert_eq!(code, 0, "{err}");
    assert_eq!(out, "{\"type\":\"FeatureCollection\",\"features\":[]}\n");
    assert!(
        err.contains("1 zone has no area (every vertex on a pole) and was left out"),
        "{err}"
    );
    let (_, out, err) = run(
        &["igeo7", "zone", "90,0", "19", "-f", "geojson", "-centroids"],
        "",
    );
    assert_eq!(out.matches("\"type\":\"Point\"").count(), 1, "{out}");
    assert!(err.is_empty(), "{err}");
    let (code, out, _) = run(&["igeo7", "zone", "90,0", "19"], "");
    assert_eq!(code, 0);
    assert!(out.contains("IGEO7 zone"), "{out}");
}

#[test]
fn text_answers_stream() {
    /// Takes what it is given until it holds more than a few kilobytes, then closes as a pipe does.
    struct Reader(usize);
    impl std::io::Write for Reader {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.0 > 4096 {
                return Err(std::io::ErrorKind::BrokenPipe.into());
            }
            self.0 += bytes.len();
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    // The whole answer is some megabytes; a run that streams stops after the first few kilobytes.
    for (zone, stdin) in [("006415654636", ""), ("-", "006415654636\n")] {
        let args: Vec<String> = ["igeo7", "disk", zone, "200", "-f", "text"]
            .iter()
            .map(|a| a.to_string())
            .collect();
        let (mut reader, mut err) = (Reader(0), Vec::new());
        let mut input = std::io::Cursor::new(stdin.as_bytes().to_vec());
        assert_eq!(
            rs4dggs_cli::run(&args, &mut input, &mut reader, &mut err),
            0
        );
        assert!(err.is_empty());
        assert!(
            reader.0 > 0 && reader.0 < 200_000,
            "{zone}: received {}",
            reader.0
        );
    }
}

#[test]
fn the_parent_that_is_a_centroid_child_is_marked_on_both_apertures() {
    let (code, out, err) = run(&["isea3h", "zone", "0.5,0.25", "7"], "");
    assert_eq!(code, 0, "{err}");
    assert!(
        out.contains("  parents     1: D4-FC-A (centroid parent)\n"),
        "{out}"
    );
    // On aperture 7 as well: of two parents the one that is itself a centroid child, in
    // whichever place the list has it; a lone parent that is one; and none where no parent is.
    for (zone, parents) in [
        (
            "00641565463601",
            "2: 0064156546360 (centroid parent), 0064156546361",
        ),
        (
            "00641565463652",
            "2: 0064156546365, 0064156546360 (centroid parent)",
        ),
        ("00641565463600", "1: 0064156546360 (centroid parent)"),
        ("0064156546360", "1: 006415654636"),
        ("006415654636", "2: 00641565463, 00641565462"),
    ] {
        let (code, out, err) = run(&["igeo7", "info", zone], "");
        assert_eq!(code, 0, "{zone}: {err}");
        assert!(
            out.contains(&format!("  parents     {parents}\n")),
            "{zone}: {out}"
        );
    }
    // The CSV names the parents, in the same order, without the mark.
    let (_, out, _) = run(&["igeo7", "info", "00641565463601", "-f", "csv"], "");
    assert_eq!(
        out.lines().nth(1).unwrap().split(',').next_back(),
        Some("0064156546360 0064156546361")
    );
}

#[test]
fn examples_are_zones_of_the_grid_asked_about() {
    let (code, out, _) = run(&["isea3h", "help", "rel"], "");
    assert_eq!(code, 0);
    assert!(out.contains("rs4dggs isea3h rel C2-23-C C4-6-B"), "{out}");
    assert_eq!(run(&["isea3h", "rel", "C2-23-C", "C4-6-B"], "").0, 0);
    let (code, _, err) = run(&["isea3h", "geom"], "");
    assert_eq!(code, 2);
    assert!(err.contains("rs4dggs isea3h geom C2-23-C"), "{err}");
    let (_, _, err) = run(&["igeo7", "sub"], "");
    assert!(err.contains("rs4dggs igeo7 sub 0064156 -depth 2"), "{err}");
    let (_, _, err) = run(&["rtea7h", "index", "0064156"], "");
    assert!(
        err.contains("rs4dggs rtea7h index 0064156 006415600"),
        "{err}"
    );
}

#[test]
fn the_examples_of_sub_and_index_run_on_every_grid() {
    for (grids, zone, depth, index, sub) in [
        (
            ["igeo7", "ivea7h", "rtea7h"],
            "0064156",
            "2",
            "27",
            "006415600",
        ),
        (["isea3h", "ivea3h", "rtea3h"], "A4-0-A", "3", "8", "B2-5-C"),
    ] {
        for g in grids {
            let (code, out, _) = run(&[g, "help", "sub"], "");
            assert_eq!(code, 0);
            for example in [
                format!("  rs4dggs {g} sub {zone} -depth {depth}\n"),
                format!("  rs4dggs {g} sub {zone} {index} -depth {depth}\n"),
            ] {
                assert!(out.contains(&example), "{g}: {out}");
            }
            let (_, out, _) = run(&[g, "help", "index"], "");
            assert!(
                out.contains(&format!("  rs4dggs {g} index {zone} {sub}\n")),
                "{g}: {out}"
            );
            // What the examples answer: a list, the sub-zone at the index, and that index back.
            let (code, out, err) = run(&[g, "sub", zone, "-depth", depth], "");
            assert_eq!(code, 0, "{g}: {err}");
            assert!(out.contains(&format!("{index:>7}  {sub}\n")), "{g}: {out}");
            let (code, out, err) = run(&[g, "sub", zone, index, "-depth", depth], "");
            assert_eq!(code, 0, "{g}: {err}");
            assert_eq!(
                out,
                format!("sub-zone {index} of {zone} at depth {depth}: {sub}\n")
            );
            let (code, out, err) = run(&[g, "index", zone, sub], "");
            assert_eq!(code, 0, "{g}: {err}");
            assert_eq!(
                out,
                format!("{sub} is sub-zone {index} of {zone}, at depth {depth}\n")
            );
        }
    }
}

#[test]
fn no_help_confines_the_hierarchy_or_the_sub_zones_to_one_aperture() {
    let (_, general, _) = run(&["help"], "");
    for line in general.lines() {
        let command = line.starts_with("  sub ") || line.starts_with("  index ");
        assert!(!(command && line.contains("aperture")), "{line}");
    }
    assert!(general.contains("\n  sub <zone> [index] "), "{general}");
    for g in ["igeo7", "ivea7h", "rtea7h", "isea3h", "ivea3h", "rtea3h"] {
        for command in ["sub", "index", "info"] {
            let (code, out, _) = run(&[g, "help", command], "");
            assert_eq!(code, 0, "{g} {command}");
            let lower = out.to_lowercase();
            for gone in ["aperture 3 only", "congruent", "defined yet", "isea3h sub"] {
                assert!(
                    g == "isea3h" && gone == "isea3h sub" || !lower.contains(gone),
                    "{g} help {command} says {gone:?}: {out}"
                );
            }
        }
    }
    // The card's help says what the hierarchy of that grid's aperture is.
    let (_, out, _) = run(&["igeo7", "help", "info"], "");
    assert!(
        out.contains("one parent or two") && out.contains("thirteen children"),
        "{out}"
    );
    let (_, out, _) = run(&["isea3h", "help", "info"], "");
    assert!(
        out.contains("one parent or three") && out.contains("seven children"),
        "{out}"
    );
}

#[test]
fn a_position_of_an_order_that_holds_no_zone_is_printed_as_no_cell() {
    // Along one of the two broken seams of the aperture-7 grids, 99 of the 17,053 positions of
    // this order hold the null zone, the first of them position 204.
    let (code, out, err) = run(&["igeo7", "sub", "010004000400", "-depth", "5"], "");
    assert_eq!(code, 0, "{err}");
    assert!(
        out.starts_with("IGEO7 zone 010004000400: 17053 sub-zones at depth 5\n"),
        "{}",
        &out[..200]
    );
    let entries: Vec<&str> = out.lines().skip(1).collect();
    assert_eq!(entries.len(), 17053);
    assert_eq!(entries[203], "    203  01000400040032330");
    assert_eq!(entries[204], "    204  (no cell)");
    assert_eq!(
        entries
            .iter()
            .filter(|l| l.ends_with("  (no cell)"))
            .count(),
        99
    );
    assert!(
        !out.contains("null"),
        "the library's own text for the null zone is printed"
    );
    let (code, out, err) = run(&["igeo7", "sub", "010004000400", "204", "-depth", "5"], "");
    assert_eq!(code, 0, "{err}");
    assert_eq!(out, "sub-zone 204 of 010004000400 at depth 5: (no cell)\n");
    // An order is written as text alone, so no other format has such a position to write.
    for format in ["csv", "geojson"] {
        let args = ["igeo7", "sub", "010004000400", "-depth", "5", "-f", format];
        let (code, out, err) = run(&args, "");
        assert_eq!((code, out.as_str()), (2, ""), "{format}: {err}");
        assert!(err.contains("sub writes text"), "{format}: {err}");
    }
}

#[test]
fn rel_follows_the_hierarchy_of_the_grid_on_aperture_7() {
    // A zone under two parents is an immediate child of each, and a sub-zone of each.
    for (parent, index) in [("00641565463", 10), ("00641565462", 1)] {
        let (code, out, err) = run(&["igeo7", "rel", parent, "006415654636"], "");
        assert_eq!(code, 0, "{err}");
        assert_eq!(
            out,
            format!(
                "{parent} is coarser than 006415654636 by 1 resolution\n\
                 {parent} and 006415654636 are not neighbours\n\
                 {parent} is an immediate parent of 006415654636\n\
                 {parent} is an ancestor of 006415654636\n\
                 006415654636 is sub-zone {index} of {parent}, at depth 1\n"
            )
        );
    }
    let (_, out, _) = run(&["igeo7", "rel", "006415654636", "00641565462"], "");
    assert!(
        out.contains("006415654636 is an immediate child of 00641565462\n")
            && out.contains("006415654636 is a descendant of 00641565462\n")
            && out.contains("006415654636 is sub-zone 1 of 00641565462, at depth 1\n"),
        "{out}"
    );
    let (_, out, _) = run(&["igeo7", "rel", "00641565463601", "00641565463602"], "");
    assert!(
        out.contains("00641565463601 and 00641565463602 are siblings\n"),
        "{out}"
    );
    // On a broken seam the engine's walk gives the index of another zone, and the library
    // finds the sub-zone in the order itself: `index` and `rel` both state its place.
    let (zone, sub) = ("00000000000000000", "000000000000000001");
    let (_, out, _) = run(&["igeo7", "sub", zone], "");
    assert!(out.contains(&format!("      3  {sub}\n")), "{out}");
    let (code, out, _) = run(&["igeo7", "index", zone, sub], "");
    assert_eq!(
        (code, out),
        (0, format!("{sub} is sub-zone 3 of {zone}, at depth 1\n"))
    );
    let (_, out, _) = run(&["igeo7", "rel", zone, sub], "");
    assert!(
        out.contains(&format!("{zone} is an immediate parent of {sub}\n"))
            && out.contains(&format!("{sub} is sub-zone 3 of {zone}, at depth 1\n")),
        "{out}"
    );
    // Where the order is longer than the library lists, the index is left unsaid as well.
    let (code, out, err) = run(&["igeo7", "rel", "00", "0000000000000"], "");
    assert_eq!(code, 0, "{err}");
    assert!(
        out.contains("00 is an ancestor of 0000000000000\n") && !out.contains("sub-zone"),
        "{out}"
    );
}

#[test]
fn a_point_with_a_space_after_the_comma_is_read() {
    let (code, out, _) = run(&["igeo7", "zone", "38.72,", "-9.14", "10"], "");
    assert_eq!(code, 0);
    assert!(out.contains("006415654636"), "{out}");
    let (code, out, _) = run(&["igeo7", "zone", "38.72,", "-9.14"], "");
    assert_eq!(code, 0);
    assert!(out.contains("006415654636"), "{out}");
    // A complete point followed by a comma is one word, not the first half of a pasted point.
    let (code, out, err) = run(&["igeo7", "zone", "38.72,-9.14,", "10"], "");
    assert_eq!(code, 0, "{err}");
    assert!(out.contains("IGEO7 zone 006415654636"), "{out}");
}

#[test]
fn negative_zero_is_zero_only_under_precision() {
    let (_, out, _) = run(
        &[
            "igeo7",
            "zone",
            "-0.0001,-0.0001",
            "3",
            "-f",
            "csv",
            "-precision",
            "3",
        ],
        "",
    );
    assert!(out.contains("\n0.000,0.000,"), "{out}");
    let (_, out, _) = run(&["igeo7", "zone", "-0,-0", "3", "-f", "csv"], "");
    assert!(out.contains("\n-0,-0,"), "{out}");
}

#[test]
fn friendlier_refusals() {
    let (code, _, err) = run(&["igeo7", "zone", "1,2", "300"], "");
    assert_eq!(code, 2);
    assert!(err.contains("0 to 255"), "{err}");
    let (_, _, err) = run(&["isea3h", "sub", "C2-23-C", "-depth", "300"], "");
    assert!(err.contains("0 to 255"), "{err}");
    let (_, _, err) = run(&["igeo7", "zone", "1,2", "3", "4"], "");
    assert!(
        err.contains("unwanted: 4") && err.contains("for example"),
        "{err}"
    );
    for topic in ["igeo7", "grids", "help"] {
        assert_eq!(run(&["help", topic], "").0, 0, "{topic}");
    }
    let (_, _, err) = run(&["igeo7", "info", "-f", "csv"], "");
    assert!(err.contains("info <zone>"), "{err}");
    let (code, _, err) = run(&["igeo7", "info", "006415654636", "-o", "-f", "csv"], "");
    assert_eq!(code, 2);
    assert!(err.contains("not the option -f"), "{err}");
    let (code, _, _) = run(
        &["isea3h", "index", "C2-23-C", "C2-23-C0", "-depth", "1"],
        "",
    );
    assert_eq!(code, 2);
    let (code, _, err) = run(&["igeo7", "info", "006415654636", "-centroids"], "");
    assert_eq!(code, 2);
    assert!(err.contains("-centroids"), "{err}");
    let (_, out, _) = run(&["isea3h", "sub", "C2-23-C", "-depth", "0"], "");
    assert!(out.contains("1 sub-zone at depth 0"), "{out}");
}

/// The features of a GeoJSON answer, each as its text (the members stay in the order written).
fn features(json: &str) -> Vec<&str> {
    json.split(r#"{"type":"Feature""#).skip(1).collect()
}

#[test]
fn neighbour_features_name_the_zone_asked_about() {
    let zones = "006415654636\n006415654634\n";
    let (code, csv, _) = run(&["igeo7", "neighbours", "-", "-f", "csv"], zones);
    assert_eq!(code, 0);
    let (code, out, _) = run(&["igeo7", "neighbours", "-", "-f", "geojson"], zones);
    assert_eq!(code, 0);
    for z in ["006415654636", "006415654634"] {
        let wanted = csv
            .lines()
            .filter(|l| l.starts_with(&format!("{z},")))
            .count();
        let mark = format!(r#""neighbour_of":"{z}""#);
        let got = features(&out).iter().filter(|f| f.contains(&mark)).count();
        assert!(wanted > 0);
        assert_eq!(got, wanted, "{z}");
    }
    assert_eq!(
        features(&out)
            .iter()
            .filter(|f| f.contains("neighbour_of"))
            .count(),
        features(&out).len()
    );
    let (_, one, _) = run(
        &["igeo7", "neighbours", "006415654636", "-f", "geojson"],
        "",
    );
    assert!(
        one.contains(r#""resolution":10,"neighbour_of":"006415654636""#),
        "{one}"
    );
}

#[test]
fn disk_features_name_their_centre() {
    let (code, out, _) = run(
        &["igeo7", "disk", "-", "1", "-f", "geojson"],
        "006415654636\n006415654634\n",
    );
    assert_eq!(code, 0);
    let all = features(&out);
    assert_eq!(all.len(), 14);
    for (i, f) in all.iter().enumerate() {
        let centre = ["006415654636", "006415654634"][i / 7];
        assert!(f.contains(&format!(r#""centre":"{centre}""#)), "{i}: {f}");
        assert!(f.contains(r#""ring":"#), "{i}: {f}");
    }
}

#[test]
fn zone_features_carry_the_point_they_answer() {
    let (code, out, _) = run(
        &[
            "igeo7",
            "zone",
            "-",
            "5",
            "-f",
            "geojson",
            "-precision",
            "2",
        ],
        "38.72,-9.14\n0 0\n",
    );
    assert_eq!(code, 0);
    let all = features(&out);
    assert_eq!(all.len(), 2);
    assert!(all[0].contains(r#""lat":38.72,"lon":-9.14"#), "{}", all[0]);
    assert!(all[1].contains(r#""lat":0.00,"lon":0.00"#), "{}", all[1]);
    // Without a resolution: every feature of a point carries it.
    let (_, out, _) = run(&["igeo7", "zone", "38.72,-9.14", "-f", "geojson"], "");
    assert_eq!(features(&out).len(), 20);
    assert!(
        features(&out)
            .iter()
            .all(|f| f.contains(r#""lat":38.72,"lon":-9.14"#))
    );
}

#[test]
fn options_that_print_no_coordinate_are_refused() {
    let cases: [&[&str]; 9] = [
        &["grids", "-centroids"],
        &["grids", "-precision", "3"],
        &["igeo7", "info", "-precision", "3"],
        &["igeo7", "-centroids"],
        &["isea3h", "rel", "C2-23-C", "C4-6-B", "-precision", "2"],
        &["isea3h", "sub", "C2-23-C", "-centroids"],
        &["isea3h", "index", "C2-23-C", "C2-23-C0", "-precision", "2"],
        &["isea3h", "rel", "C2-23-C", "C4-6-B", "-centroids"],
        &["isea3h", "index", "C2-23-C", "C2-23-C0", "-centroids"],
    ];
    for args in cases {
        let (code, out, err) = run(args, "");
        assert_eq!(code, 2, "{args:?}: {err}");
        assert!(out.is_empty(), "{args:?}: {out}");
        assert!(err.contains("no coordinate"), "{args:?}: {err}");
    }
    // Where a coordinate is printed, both are still accepted.
    let (code, _, err) = run(&["igeo7", "info", "006415654636", "-precision", "3"], "");
    assert_eq!(code, 0, "{err}");
    let (_, _, err) = run(&["isea3h", "sub", "C2-23-C", "-centroids"], "");
    assert!(!err.contains("add -f geojson"), "{err}");
}

#[test]
fn help_names_a_grid_then_a_command() {
    let (code, out, _) = run(&["help", "isea3h", "sub"], "");
    assert_eq!(code, 0);
    assert!(out.contains("rs4dggs isea3h sub A4-0-A -depth 3"), "{out}");
    let (code, out, _) = run(&["help", "igeo7", "zone"], "");
    assert_eq!(code, 0);
    assert!(out.contains("rs4dggs igeo7 zone 38.72,-9.14 10"), "{out}");
    let (code, out, _) = run(&["help", "igeo7", "neighbors"], "");
    assert_eq!(code, 0);
    assert!(out.contains("igeo7 neighbours"), "{out}");
    let (code, _, err) = run(&["help", "igeo7", "nonsense"], "");
    assert_eq!(code, 2);
    assert!(err.contains("unknown command"), "{err}");
    // `help igeo7` alone is still that grid's info.
    let (_, out, _) = run(&["help", "igeo7"], "");
    assert!(out.starts_with("IGEO7\n"), "{out}");
}

#[test]
fn the_helps_say_that_a_first_line_may_be_a_header() {
    let (_, general, _) = run(&["help"], "");
    let (_, zone, _) = run(&["igeo7", "help", "zone"], "");
    let (_, info, _) = run(&["igeo7", "help", "info"], "");
    for (what, text) in [("help", general), ("zone", zone), ("info", info)] {
        assert!(text.contains("taken as a header"), "{what}: {text}");
    }
}

#[test]
fn no_line_of_the_tool_exceeds_a_hundred_columns() {
    let dir = format!("{}/src", env!("CARGO_MANIFEST_DIR"));
    let mut checked = 0;
    for entry in std::fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_some_and(|e| e == "rs") {
            let source = std::fs::read_to_string(&path).unwrap();
            for (n, line) in source.lines().enumerate() {
                assert!(
                    line.chars().count() <= 100,
                    "{}:{}: {line}",
                    path.display(),
                    n + 1
                );
            }
            checked += 1;
        }
    }
    assert!(checked >= 8, "only {checked} source files found in {dir}");
}

#[test]
fn help_for_a_grid_is_never_refused() {
    let (_, plain, _) = run(&["help", "igeo7"], "");
    assert!(plain.starts_with("IGEO7\n"), "{plain}");
    for extra in [&["-precision", "3"][..], &["-centroids"], &["-f", "csv"]] {
        let mut args = vec!["help", "igeo7"];
        args.extend_from_slice(extra);
        let (code, out, err) = run(&args, "");
        assert_eq!(code, 0, "{args:?}: {err}");
        assert_eq!(out, plain, "{args:?}");
    }
    // Not a request for help: still refused.
    assert_eq!(run(&["igeo7", "info", "-precision", "3"], "").0, 2);
}

/// The tool as a process of its own, with the limit on a list of sub-zones lowered to 20 by
/// its environment: the variable is read once in a process, so that it cannot be set for one
/// test of a process that runs many.
fn run_with_a_limit_of_20(args: &[&str]) -> (Option<i32>, String, String) {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_rs4dggs"))
        .args(args)
        .env("RS4DGGS_MAX_MATERIALISED_SUB_ZONES", "20")
        .output()
        .unwrap();
    (
        output.status.code(),
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8(output.stderr).unwrap(),
    )
}

#[test]
fn the_environment_lowers_the_limit_on_a_list_of_sub_zones() {
    // On an aperture-7 grid: 13 sub-zones are within the limit, and 55 are refused with the
    // limit in force; so is the index of a zone in an order of 55, while one in an order of
    // 13 is answered, and so is one sub-zone by its index at any depth.
    let (code, out, _) = run_with_a_limit_of_20(&["igeo7", "sub", "0064156"]);
    assert_eq!(code, Some(0));
    assert!(out.contains("13 sub-zones at depth 1"), "{out}");
    let (code, _, err) = run_with_a_limit_of_20(&["igeo7", "sub", "0064156", "-depth", "2"]);
    assert_eq!(code, Some(1));
    assert!(
        err.contains("55 sub-zones are more than the limit of 20"),
        "{err}"
    );
    let (code, out, _) = run_with_a_limit_of_20(&["igeo7", "index", "0064156", "00641565"]);
    assert_eq!(code, Some(0));
    assert!(out.contains("is sub-zone 5 of 0064156"), "{out}");
    let (code, _, err) = run_with_a_limit_of_20(&["igeo7", "index", "0064156", "006415600"]);
    assert_eq!(code, Some(1));
    assert!(
        err.contains("55 sub-zones are more than the limit of 20"),
        "{err}"
    );
    let (code, out, _) = run_with_a_limit_of_20(&["igeo7", "sub", "0064156", "27", "-depth", "2"]);
    assert_eq!(code, Some(0));
    assert!(out.contains("006415600"), "{out}");
    // An order longer than the ceiling itself is refused with the limit in force too, for
    // the list and for an index.
    for args in [
        &["igeo7", "sub", "0064156", "-depth", "8"][..],
        &["igeo7", "index", "0064156", "006415600000000"],
    ] {
        let (code, _, err) = run_with_a_limit_of_20(args);
        assert_eq!(code, Some(1), "{args:?}");
        assert!(
            err.contains("5767201 sub-zones are more than the limit of 20"),
            "{args:?}: {err}"
        );
    }
    // On an aperture-3 grid: 6 are within it, and 31 are refused.
    let (code, out, _) = run_with_a_limit_of_20(&["isea3h", "sub", "A4-0-A"]);
    assert_eq!(code, Some(0));
    assert!(out.contains("6 sub-zones at depth 1"), "{out}");
    let (code, _, err) = run_with_a_limit_of_20(&["isea3h", "sub", "A4-0-A", "-depth", "3"]);
    assert_eq!(code, Some(1));
    assert!(
        err.contains("31 sub-zones are more than the limit of 20"),
        "{err}"
    );
    let (code, _, err) = run_with_a_limit_of_20(&["isea3h", "index", "A4-0-A", "B2-5-C"]);
    assert_eq!(code, Some(1));
    assert!(
        err.contains("31 sub-zones are more than the limit of 20"),
        "{err}"
    );
}
