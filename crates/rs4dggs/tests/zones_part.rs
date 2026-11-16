//! A level's walk cut into contiguous parts (`zones_part`), for a caller who walks one level on
//! several threads: the parts, taken in order, are the walk of the level, zone for zone.
use rs4dggs::{AnyGrid, Error, ZoneId};

const GRIDS: [&str; 6] = ["ISEA3H", "IVEA3H", "RTEA3H", "IGEO7", "IVEA7H", "RTEA7H"];

fn grid(name: &str) -> AnyGrid {
    rs4dggs::get_grid(name).unwrap()
}

fn part(g: AnyGrid, level: u8, k: u64, n: u64) -> Vec<ZoneId> {
    g.zones_part(level, k, n).unwrap().collect()
}

#[test]
fn the_parts_in_order_are_the_walk() {
    for name in GRIDS {
        let g = grid(name);
        for level in [0, 1, 3] {
            let whole: Vec<ZoneId> = g.zones(level).unwrap().collect();
            // 64 parts at level 0 is more parts than rows: some are empty, none is lost.
            for n in [1, 2, 3, 7, 13, 64] {
                let cut: Vec<ZoneId> = (0..n).flat_map(|k| part(g, level, k, n)).collect();
                assert_eq!(cut, whole, "{name} level {level} in {n} parts");
            }
        }
    }
}

#[test]
fn the_parts_hold_about_as_many_zones() {
    for (name, level) in [("IGEO7", 6), ("ISEA3H", 9)] {
        let g = grid(name);
        let sizes: Vec<usize> = (0..8)
            .map(|k| g.zones_part(level, k, 8).unwrap().count())
            .collect();
        assert_eq!(
            sizes.iter().sum::<usize>() as u64,
            g.count_zones(level).unwrap()
        );
        let (min, max) = (*sizes.iter().min().unwrap(), *sizes.iter().max().unwrap());
        assert!(
            max as f64 <= 1.25 * min as f64,
            "{name} level {level}: {sizes:?}"
        );
    }
}

#[test]
fn a_part_entered_after_one_of_its_zones_ends_with_it() {
    let g = grid("IGEO7");
    let third = part(g, 4, 1, 3);
    assert!(third.len() > 20);
    let rest: Vec<ZoneId> = g
        .zones_part(4, 1, 3)
        .unwrap()
        .after(third[9])
        .unwrap()
        .collect();
    assert_eq!(rest, third[10..]);
}

#[test]
fn a_part_out_of_range_is_refused() {
    let g = grid("IGEO7");
    assert_eq!(
        g.zones_part(2, 3, 3).err(),
        Some(Error::PartOutOfRange { part: 3, parts: 3 })
    );
    assert_eq!(
        g.zones_part(2, 0, 0).err(),
        Some(Error::PartOutOfRange { part: 0, parts: 0 })
    );
    assert!(matches!(
        g.zones_part(g.max_box_level() + 1, 0, 1),
        Err(Error::ResolutionOutOfRange { .. })
    ));
}

#[test]
fn a_typed_grid_cuts_the_same_parts() {
    let typed = rs4dggs::igeo7();
    let any = grid("IGEO7");
    for k in 0..4 {
        let a: Vec<ZoneId> = typed.zones_part(3, k, 4).unwrap().collect();
        assert_eq!(a, part(any, 3, k, 4));
    }
}
