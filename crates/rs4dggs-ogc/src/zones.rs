//! The list of zones that answers a zone query: as JSON, and as 64-bit integers.

use std::fmt;
use std::io::Write;

use rs4dggs::{AnyGrid, ZoneId};

use crate::link::{self, Link, LinkTemplate};
use crate::{Error, Result, json};

/// The writer of a list of zones as JSON, the document a zone query answers with
/// `application/json`: `{"zones":[...],"returnedAreaMetersSquare":...,"links":[...],
/// "linkTemplates":[...]}`.
///
/// The zones are given one at a time, in the order they are to appear, and each is written as
/// the grid's text identifier as soon as it is given; [`ZoneListJson::finish`] writes the rest
/// of the document. No white space is written, nor a final line feed.
///
/// The writer does not judge the links, which are the service's. The standard requires two in
/// this document: one to the description of the grid, with the relation type [`rel::DGGRS`],
/// and one to its definition, with [`rel::DGGRS_DEFINITION`].
///
/// A failure of the sink in the middle of the list leaves the document truncated where it
/// arose.
///
/// [`rel::DGGRS`]: crate::rel::DGGRS
/// [`rel::DGGRS_DEFINITION`]: crate::rel::DGGRS_DEFINITION
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{Link, ZoneListJson, rel};
///
/// let grid = rs4dggs::get_grid("ISEA3H")?;
/// let mut list = ZoneListJson::new(Vec::new(), grid);
/// list.zone(grid.zone_from_text("A6-0-C")?)?;
/// let links = [
///     Link { rel: rel::DGGRS, href: "/dggs/ISEA3H", title: None, media_type: None },
///     Link { rel: rel::DGGRS_DEFINITION, href: "/dggrs/ISEA3H", title: None, media_type: None },
/// ];
/// assert_eq!(list.finish(None, &links, &[])?, 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct ZoneListJson<W: Write> {
    out: W,
    grid: AnyGrid,
    zones: u64,
}

impl<W: Write> ZoneListJson<W> {
    /// A list of the zones of `grid`, to be written to `out`. Nothing is written yet.
    pub fn new(out: W, grid: AnyGrid) -> Self {
        Self {
            out,
            grid,
            zones: 0,
        }
    }

    /// Writes one zone, as its text identifier.
    ///
    /// The null zone is refused with [`Error::NullZone`] before anything is written: it names
    /// no zone, and a list that may hold it is filtered by the caller.
    ///
    /// The zone must be one the grid defines, as the library gives them (`zone_from_text`,
    /// `zones`, `zones_in_box`, `compact_zones`, `sub_zones`). Nothing else is checked: an
    /// identifier the grid does not define is written as the library's `text_id` renders it,
    /// which may be the text of another zone.
    pub fn zone(&mut self, zone: ZoneId) -> Result<()> {
        if zone == ZoneId::NULL {
            return Err(Error::NullZone);
        }
        self.out.write_all(if self.zones == 0 {
            b"{\"zones\":["
        } else {
            b","
        })?;
        json::string(&mut self.out, &self.grid.text_id(zone))?;
        self.zones += 1;
        Ok(())
    }

    /// Ends the document and returns the number of zones written.
    ///
    /// `returned_area` is written as `returnedAreaMetersSquare` when given and finite; a value
    /// that is not finite is left out, as if it were not given, since the member is optional
    /// and must be a number where present. The standard asks that it be the area of the list
    /// before any compaction, so that zones a compaction would make overlap count once: the sum
    /// of the library's `area` over the zones of the uncompacted list. `links` are always
    /// written, as an empty array when there are none; `templates` only when there is at least
    /// one. The sink is flushed.
    pub fn finish(
        mut self,
        returned_area: Option<f64>,
        links: &[Link<'_>],
        templates: &[LinkTemplate<'_>],
    ) -> Result<u64> {
        let out = &mut self.out;
        out.write_all(if self.zones == 0 {
            b"{\"zones\":[]"
        } else {
            b"]"
        })?;
        if let Some(area) = returned_area.filter(|a| a.is_finite()) {
            out.write_all(b",\"returnedAreaMetersSquare\":")?;
            json::f64(out, area)?;
        }
        out.write_all(b",\"links\":")?;
        link::links(out, links)?;
        if !templates.is_empty() {
            out.write_all(b",\"linkTemplates\":")?;
            link::templates(out, templates)?;
        }
        out.write_all(b"}")?;
        out.flush()?;
        Ok(self.zones)
    }
}

impl<W: Write> fmt::Debug for ZoneListJson<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ZoneListJson")
            .field("grid", &self.grid.name())
            .field("zones", &self.zones)
            .finish_non_exhaustive()
    }
}

/// Writes a list of zones as 64-bit integers, the document a zone query answers with
/// `application/x-binary`: the number of zones, then each zone's identifier (`ZoneId.0`), every
/// one of them eight bytes, least significant byte first. Nothing else is written: no header,
/// no padding, no trailer. The sink is flushed.
///
/// The identifiers are unsigned. Those of the zones under the zones 08 to 11 of level 0 of the
/// three grids of aperture 7 have their most significant bit set, so that a client that reads
/// the eight bytes as a signed integer sees a negative number for them.
///
/// A slice that holds the null zone is refused with [`Error::NullZone`] before anything is
/// written. The zones must be ones the grid defines, as the library gives them; nothing else is
/// checked, and an identifier the grid does not define is written as it is. A failure of the
/// sink leaves the document truncated where it arose.
pub fn write_zone_list_uint64<W: Write>(mut out: W, zones: &[ZoneId]) -> Result<()> {
    if zones.contains(&ZoneId::NULL) {
        return Err(Error::NullZone);
    }
    out.write_all(&(zones.len() as u64).to_le_bytes())?;
    for zone in zones {
        out.write_all(&zone.0.to_le_bytes())?;
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    fn isea3h() -> AnyGrid {
        rs4dggs::get_grid("ISEA3H").unwrap()
    }

    fn id(grid: AnyGrid, text: &str) -> ZoneId {
        grid.zone_from_text(text).unwrap()
    }

    /// A sink that accepts ten bytes and refuses everything after them.
    struct TenBytes(Vec<u8>);

    impl Write for TenBytes {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            let room = 10 - self.0.len();
            if room == 0 {
                return Err(io::Error::other("full"));
            }
            let n = room.min(buf.len());
            self.0.extend_from_slice(&buf[..n]);
            Ok(n)
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn an_empty_list_has_its_zones_and_its_links() {
        let mut out = Vec::new();
        let list = ZoneListJson::new(&mut out, isea3h());
        assert_eq!(list.finish(None, &[], &[]).unwrap(), 0);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "{\"zones\":[],\"links\":[]}"
        );
    }

    #[test]
    fn an_area_that_is_not_finite_is_left_out() {
        for area in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut out = Vec::new();
            let list = ZoneListJson::new(&mut out, isea3h());
            list.finish(Some(area), &[], &[]).unwrap();
            assert_eq!(
                String::from_utf8(out).unwrap(),
                "{\"zones\":[],\"links\":[]}",
                "{area}"
            );
        }
    }

    #[test]
    fn the_area_links_and_templates_follow_the_zones_in_order() {
        let grid = isea3h();
        let mut out = Vec::new();
        let mut list = ZoneListJson::new(&mut out, grid);
        list.zone(id(grid, "A6-0-C")).unwrap();
        list.zone(id(grid, "AA-0-B")).unwrap();
        let link = Link {
            rel: crate::rel::DGGRS,
            href: "/dggs/ISEA3H",
            title: Some("The \"ISEA3H\" grid"),
            media_type: None,
        };
        let template = LinkTemplate {
            rel: crate::rel::DGGRS_ZONE_DATA,
            uri_template: "/dggs/ISEA3H/zones/{zoneId}/data",
            title: None,
        };
        assert_eq!(list.finish(Some(0.5), &[link], &[template]).unwrap(), 2);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "{\"zones\":[\"A6-0-C\",\"AA-0-B\"],\"returnedAreaMetersSquare\":0.5,\
             \"links\":[{\"rel\":\"https://www.opengis.net/def/rel/ogc/1.0/dggrs\",\
             \"title\":\"The \\\"ISEA3H\\\" grid\",\"href\":\"/dggs/ISEA3H\"}],\
             \"linkTemplates\":[{\"rel\":\"https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-data\",\
             \"uriTemplate\":\"/dggs/ISEA3H/zones/{zoneId}/data\"}]}"
        );
    }

    #[test]
    fn the_null_zone_is_refused_and_nothing_is_written() {
        let grid = isea3h();
        let mut out = Vec::new();
        {
            let mut list = ZoneListJson::new(&mut out, grid);
            assert!(matches!(list.zone(ZoneId::NULL), Err(Error::NullZone)));
        }
        assert!(out.is_empty());

        let mut list = ZoneListJson::new(&mut out, grid);
        list.zone(id(grid, "A6-0-C")).unwrap();
        assert!(matches!(list.zone(ZoneId::NULL), Err(Error::NullZone)));
        assert_eq!(list.finish(None, &[], &[]).unwrap(), 1);
        assert_eq!(out, b"{\"zones\":[\"A6-0-C\"],\"links\":[]}");
    }

    #[test]
    fn a_sink_that_fails_is_an_error_and_no_panic() {
        let grid = isea3h();
        let mut list = ZoneListJson::new(TenBytes(Vec::new()), grid);
        // `{"zones":[` is ten bytes: the first zone is the first write refused.
        assert!(matches!(list.zone(id(grid, "A6-0-C")), Err(Error::Io(_))));
        assert!(matches!(list.finish(None, &[], &[]), Err(Error::Io(_))));
        let empty = ZoneListJson::new(TenBytes(Vec::new()), grid);
        assert!(matches!(empty.finish(None, &[], &[]), Err(Error::Io(_))));
        let ids = [id(grid, "A6-0-C"), id(grid, "AA-0-B")];
        assert!(matches!(
            write_zone_list_uint64(TenBytes(Vec::new()), &ids),
            Err(Error::Io(_))
        ));
    }

    #[test]
    fn the_writer_shows_its_grid_and_count() {
        let list = ZoneListJson::new(Vec::new(), isea3h());
        assert_eq!(
            format!("{list:?}"),
            "ZoneListJson { grid: \"ISEA3H\", zones: 0, .. }"
        );
    }

    #[test]
    fn the_binary_list_is_a_count_and_the_identifiers_little_endian() {
        let grid = isea3h();
        let mut out = Vec::new();
        write_zone_list_uint64(&mut out, &[id(grid, "A6-0-C"), id(grid, "AA-0-B")]).unwrap();
        assert_eq!(
            out,
            [
                0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, //
                0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc0, 0x00, //
                0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x01,
            ]
        );
    }

    #[test]
    fn the_binary_list_keeps_bit_63_of_the_aperture_7_identifiers() {
        let grid = rs4dggs::get_grid("IGEO7").unwrap();
        let mut out = Vec::new();
        write_zone_list_uint64(&mut out, &[id(grid, "0064"), id(grid, "08")]).unwrap();
        assert_eq!(
            out,
            [
                0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, //
                0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f, 0x0d, //
                0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x8f,
            ]
        );
    }

    #[test]
    fn an_empty_binary_list_is_a_count_of_nought() {
        let mut out = Vec::new();
        write_zone_list_uint64(&mut out, &[]).unwrap();
        assert_eq!(out, [0u8; 8]);
    }

    #[test]
    fn a_binary_list_with_the_null_zone_is_refused_and_nothing_is_written() {
        let grid = isea3h();
        let mut out = Vec::new();
        let zones = [id(grid, "A6-0-C"), ZoneId::NULL];
        assert!(matches!(
            write_zone_list_uint64(&mut out, &zones),
            Err(Error::NullZone)
        ));
        assert!(out.is_empty());
    }
}
