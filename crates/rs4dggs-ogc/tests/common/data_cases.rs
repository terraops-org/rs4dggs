// The data of the goldens `tests/goldens/data-*.json`, which the tests of both encodings of zone
// data share: the writer of DGGS-JSON must give each golden, and the writer of DGGS-UBJSON a
// document that reads back as it. The crate that includes this file has `Dggrs`, `Property`,
// `Values` and `ZoneData` at its root.

use crate::{Dggrs, Property, Values, ZoneData};

/// Calls `f` with the name of each golden and the data it holds.
pub fn each(mut f: impl FnMut(&str, &ZoneData<'_>)) {
    let isea3h = Dggrs::from_id("ISEA3H").unwrap();
    let z7 = Dggrs::from_id("ISEA7H_Z7").unwrap();
    let c2 = isea3h.grid().zone_from_text("C2-23-C").unwrap();
    let positions = |n: i32| (0..n).map(Some).collect::<Vec<_>>();
    let seven = positions(7);
    let mut null_at_2 = seven.clone();
    null_at_2[2] = None;
    let thirteen = positions(13);

    f(
        "data-ISEA3H-C2-23-C-depth1.json",
        &ZoneData::new(
            isea3h,
            c2,
            &[1],
            &[Property::new("i", &[Values::I32(&seven)])],
        ),
    );
    f(
        "data-ISEA3H-C2-23-C-depth1.null-at-2.json",
        &ZoneData::new(
            isea3h,
            c2,
            &[1],
            &[Property::new("i", &[Values::I32(&null_at_2)])],
        ),
    );
    f(
        "data-ISEA3H-C2-23-C-depths1-2.json",
        &ZoneData::new(
            isea3h,
            c2,
            &[1, 2],
            &[Property::new(
                "i",
                &[Values::I32(&seven), Values::I32(&thirteen)],
            )],
        ),
    );
    f(
        "data-ISEA3H-C2-23-C-depth0.json",
        &ZoneData::new(
            isea3h,
            c2,
            &[0],
            &[Property::new("i", &[Values::U16(&[Some(1000)])])],
        ),
    );
    let t = [
        Some(0.5),
        Some(1.5),
        None,
        Some(0.1),
        Some(22.8),
        Some(f32::NAN),
        Some(6.0),
    ];
    f(
        "data-ISEA3H-C2-23-C-depth1.f32.json",
        &ZoneData::new(isea3h, c2, &[1], &[Property::new("t", &[Values::F32(&t)])]),
    );
    let class = [
        Some(0),
        Some(1),
        Some(2),
        Some(3),
        Some(4),
        Some(5),
        Some(255),
    ];
    let height = [
        Some(-12.5),
        Some(-0.0),
        Some(1e-7),
        Some(1e17),
        Some(1234.5678),
        Some(f64::INFINITY),
        None,
    ];
    f(
        "data-ISEA3H-C2-23-C-depth1.two-bands.json",
        &ZoneData::new(
            isea3h,
            c2,
            &[1],
            &[
                Property::new("class", &[Values::U8(&class)]),
                Property::new("height", &[Values::F64(&height)]),
            ],
        ),
    );
    let z0064 = z7.grid().zone_from_text("0064").unwrap();
    let fifty_five = positions(55);
    f(
        "data-ISEA7H_Z7-0064-depth2.json",
        &ZoneData::new(
            z7,
            z0064,
            &[2],
            &[Property::new("i", &[Values::I32(&fifty_five)])],
        ),
    );
    // At this zone, on a broken seam, two entries of the order are the null zone: no zone, and
    // so no value.
    let seam = z7.grid().zone_from_text("00055353260226021").unwrap();
    let at_seam: Vec<Option<i32>> = z7
        .grid()
        .sub_zones(seam, 1)
        .unwrap()
        .iter()
        .zip(0..)
        .map(|(&z, i)| (z != rs4dggs::ZoneId::NULL).then_some(i))
        .collect();
    f(
        "data-ISEA7H_Z7-00055353260226021-depth1.seam.json",
        &ZoneData::new(
            z7,
            seam,
            &[1],
            &[Property::new("i", &[Values::I32(&at_seam)])],
        ),
    );
}
