//! Links and link templates, as the documents of OGC API - DGGS carry them, and the relation
//! types the standard defines for them.

use std::io::{self, Write};

use crate::json;

/// A link of a document: its relation type, its target, and optionally a title and the media
/// type of the target.
///
/// The writers take the links as the service gives them and judge none of them.
///
/// Built by [`Link::new`], with [`Link::with_title`] and [`Link::with_media_type`] for the
/// members that may be absent; the type may gain fields in a later release, so that it is not
/// built from its fields outside this crate.
///
/// ```
/// use rs4dggs_ogc::{Link, rel};
///
/// let link = Link::new(rel::DGGRS_DEFINITION, "/dggrs/ISEA3H").with_media_type("application/json");
/// assert_eq!((link.title, link.media_type), (None, Some("application/json")));
/// ```
///
/// Built from its fields outside this crate, it does not compile:
///
/// ```compile_fail,E0639
/// let link = rs4dggs_ogc::Link { rel: "self", href: "/", title: None, media_type: None };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct Link<'a> {
    /// The relation type, in full: one of [`rel`] or any other.
    pub rel: &'a str,
    /// The target.
    pub href: &'a str,
    /// A title for the target, written when present.
    pub title: Option<&'a str>,
    /// The media type of the target, written as the member `type` when present.
    pub media_type: Option<&'a str>,
}

/// A link template: a link whose target holds variables, such as `{zoneId}`.
///
/// Built by [`LinkTemplate::new`], with [`LinkTemplate::with_title`] for a title; the type may
/// gain fields in a later release, so that it is not built from its fields outside this crate.
///
/// ```
/// use rs4dggs_ogc::{LinkTemplate, rel};
///
/// let template = LinkTemplate::new(rel::DGGRS_ZONE_DATA, "/dggs/ISEA3H/zones/{zoneId}/data")
///     .with_title("Data of a zone");
/// assert_eq!(template.title, Some("Data of a zone"));
/// ```
///
/// Built from its fields outside this crate, it does not compile:
///
/// ```compile_fail,E0639
/// let template = rs4dggs_ogc::LinkTemplate { rel: "self", uri_template: "/{zoneId}", title: None };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct LinkTemplate<'a> {
    /// The relation type, in full.
    pub rel: &'a str,
    /// The target, with its variables between braces.
    pub uri_template: &'a str,
    /// A title for the target, written when present.
    pub title: Option<&'a str>,
}

impl<'a> Link<'a> {
    /// A link of relation type `rel` to `href`, without a title or a media type.
    pub const fn new(rel: &'a str, href: &'a str) -> Self {
        Self {
            rel,
            href,
            title: None,
            media_type: None,
        }
    }

    /// The same link with the title `title`.
    pub const fn with_title(self, title: &'a str) -> Self {
        Self {
            title: Some(title),
            ..self
        }
    }

    /// The same link with the media type `media_type`.
    pub const fn with_media_type(self, media_type: &'a str) -> Self {
        Self {
            media_type: Some(media_type),
            ..self
        }
    }
}

impl<'a> LinkTemplate<'a> {
    /// A link template of relation type `rel` to `uri_template`, without a title.
    pub const fn new(rel: &'a str, uri_template: &'a str) -> Self {
        Self {
            rel,
            uri_template,
            title: None,
        }
    }

    /// The same template with the title `title`.
    pub const fn with_title(self, title: &'a str) -> Self {
        Self {
            title: Some(title),
            ..self
        }
    }
}

/// The relation types of OGC API - DGGS, written in full as the text of the standard's
/// requirements writes them.
pub mod rel {
    /// The description of a grid (a DGGRS) as the service offers it.
    pub const DGGRS: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs";
    /// The definition of a grid.
    pub const DGGRS_DEFINITION: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-definition";
    /// The list of the grids the service offers.
    pub const DGGRS_LIST: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-list";
    /// The information of one zone.
    pub const DGGRS_ZONE_INFO: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-info";
    /// The data of one zone.
    pub const DGGRS_ZONE_DATA: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-data";
    /// The query of the zones of a grid.
    pub const DGGRS_ZONE_QUERY: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-query";
    /// A parent of a zone.
    pub const DGGRS_ZONE_PARENT: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-parent";
    /// A child of a zone.
    pub const DGGRS_ZONE_CHILD: &str = "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-child";
    /// A neighbour of a zone.
    pub const DGGRS_ZONE_NEIGHBOR: &str =
        "https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-neighbor";
    /// The data from which the zones' values are drawn.
    pub const GEODATA: &str = "https://www.opengis.net/def/rel/ogc/1.0/geodata";
    /// The data set the resource belongs to.
    pub const DATASET: &str = "https://www.opengis.net/def/rel/ogc/1.0/dataset";
}

/// Writes `links` as a JSON array, each member in the order `rel`, `type`, `title`, `href`, the
/// absent ones left out.
pub(crate) fn links<W: Write>(out: &mut W, links: &[Link<'_>]) -> io::Result<()> {
    out.write_all(b"[")?;
    for (i, link) in links.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        out.write_all(b"{\"rel\":")?;
        json::string(out, link.rel)?;
        if let Some(media_type) = link.media_type {
            out.write_all(b",\"type\":")?;
            json::string(out, media_type)?;
        }
        if let Some(title) = link.title {
            out.write_all(b",\"title\":")?;
            json::string(out, title)?;
        }
        out.write_all(b",\"href\":")?;
        json::string(out, link.href)?;
        out.write_all(b"}")?;
    }
    out.write_all(b"]")
}

/// Writes `templates` as a JSON array, each member in the order `rel`, `title`, `uriTemplate`, an
/// absent title left out.
pub(crate) fn templates<W: Write>(out: &mut W, templates: &[LinkTemplate<'_>]) -> io::Result<()> {
    out.write_all(b"[")?;
    for (i, template) in templates.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        out.write_all(b"{\"rel\":")?;
        json::string(out, template.rel)?;
        if let Some(title) = template.title {
            out.write_all(b",\"title\":")?;
            json::string(out, title)?;
        }
        out.write_all(b",\"uriTemplate\":")?;
        json::string(out, template.uri_template)?;
        out.write_all(b"}")?;
    }
    out.write_all(b"]")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(write: impl FnOnce(&mut Vec<u8>) -> io::Result<()>) -> String {
        let mut out = Vec::new();
        write(&mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn every_relation_type_is_the_base_and_its_last_segment() {
        let base = "https://www.opengis.net/def/rel/ogc/1.0/";
        for (constant, segment) in [
            (rel::DGGRS, "dggrs"),
            (rel::DGGRS_DEFINITION, "dggrs-definition"),
            (rel::DGGRS_LIST, "dggrs-list"),
            (rel::DGGRS_ZONE_INFO, "dggrs-zone-info"),
            (rel::DGGRS_ZONE_DATA, "dggrs-zone-data"),
            (rel::DGGRS_ZONE_QUERY, "dggrs-zone-query"),
            (rel::DGGRS_ZONE_PARENT, "dggrs-zone-parent"),
            (rel::DGGRS_ZONE_CHILD, "dggrs-zone-child"),
            (rel::DGGRS_ZONE_NEIGHBOR, "dggrs-zone-neighbor"),
            (rel::GEODATA, "geodata"),
            (rel::DATASET, "dataset"),
        ] {
            assert_eq!(constant, format!("{base}{segment}"));
        }
    }

    #[test]
    fn a_link_writes_its_members_in_order_and_leaves_out_the_absent() {
        let bare = Link {
            rel: rel::DGGRS,
            href: "/dggs/ISEA3H",
            title: None,
            media_type: None,
        };
        let full = Link {
            rel: "self",
            href: "/a?f=json",
            title: Some("The \"zone\" list"),
            media_type: Some("application/json"),
        };
        assert_eq!(
            text(|o| links(o, &[bare, full])),
            "[{\"rel\":\"https://www.opengis.net/def/rel/ogc/1.0/dggrs\",\"href\":\"/dggs/ISEA3H\"},\
             {\"rel\":\"self\",\"type\":\"application/json\",\"title\":\"The \\\"zone\\\" list\",\
             \"href\":\"/a?f=json\"}]"
        );
        assert_eq!(text(|o| links(o, &[])), "[]");
    }

    #[test]
    fn a_template_writes_its_members_in_order_and_leaves_out_an_absent_title() {
        let titled = LinkTemplate {
            rel: rel::DGGRS_ZONE_DATA,
            uri_template: "/dggs/ISEA3H/zones/{zoneId}/data",
            title: Some("Data of a zone"),
        };
        let bare = LinkTemplate {
            title: None,
            ..titled
        };
        assert_eq!(
            text(|o| templates(o, &[titled, bare])),
            "[{\"rel\":\"https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-data\",\
             \"title\":\"Data of a zone\",\"uriTemplate\":\"/dggs/ISEA3H/zones/{zoneId}/data\"},\
             {\"rel\":\"https://www.opengis.net/def/rel/ogc/1.0/dggrs-zone-data\",\
             \"uriTemplate\":\"/dggs/ISEA3H/zones/{zoneId}/data\"}]"
        );
        assert_eq!(text(|o| templates(o, &[])), "[]");
    }
}
