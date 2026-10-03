//! The six grids as an OGC API names them: their identifiers, titles, URIs, coordinate reference
//! systems and definition documents.

use rs4dggs::AnyGrid;

/// One of the six grids, as a service of OGC API - DGGS names and describes it.
///
/// The identifiers are DGGAL's own, which on aperture 7 are not the library's: `ISEA7H_Z7`,
/// `IVEA7H_Z7` and `RTEA7H_Z7` are the library's `IGEO7`, `IVEA7H` and `RTEA7H`. The names
/// `ISEA7H`, `IVEA7H` and `RTEA7H` denote, in DGGAL (and, for the first two, in the OGC register),
/// grids of another indexing, and are therefore refused by [`Dggrs::from_id`]. A text received
/// from a client is resolved here and never passed to [`rs4dggs::get_grid`].
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::Dggrs;
///
/// let dggrs = Dggrs::from_id("IGEO7").unwrap();
/// assert_eq!(dggrs.id(), "ISEA7H_Z7");
/// assert_eq!(dggrs.grid().name(), "IGEO7");
/// assert!(Dggrs::from_id("IVEA7H").is_none());
/// ```
pub struct Dggrs {
    id: &'static str,
    title: &'static str,
    description: &'static str,
    uri: &'static str,
    registered: bool,
    crs: Option<&'static str>,
    grid: fn() -> AnyGrid,
    definition: &'static str,
}

/// The coordinate reference system of the two grids on the ISEA projection, as the registered
/// definitions of ISEA3H and ISEA7H give it.
const ISEA_CRS: &str = "https://www.opengis.net/def/crs/OGC/0/1534";

impl std::fmt::Debug for Dggrs {
    /// The identifier alone: the rest is the table's, and the definition is a whole document.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Dggrs")
            .field("id", &self.id)
            .finish_non_exhaustive()
    }
}

static ALL: [Dggrs; 6] = [
    Dggrs {
        id: "ISEA3H",
        title: "ISEA Aperture 3 Hexagonal",
        description: "A Discrete Global Grid Reference System based on the Icosahedral Snyder \
                      Equal Area projection, with aperture 3 hexagonal zones, using an indexing \
                      scheme based on ISEA9R.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H",
        registered: true,
        crs: Some(ISEA_CRS),
        grid: || AnyGrid::Isea3h(rs4dggs::isea3h()),
        definition: include_str!("../definitions/ISEA3H.json"),
    },
    Dggrs {
        id: "IVEA3H",
        title: "IVEA Aperture 3 Hexagonal",
        description: "A Discrete Global Grid Reference System based on the Icosahedral \
                      Vertex-oriented Equal Area projection, with aperture 3 hexagonal zones, \
                      using an indexing scheme based on IVEA9R.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/IVEA3H",
        registered: true,
        crs: None,
        grid: || AnyGrid::Ivea3h(rs4dggs::ivea3h()),
        definition: include_str!("../definitions/IVEA3H.json"),
    },
    Dggrs {
        id: "RTEA3H",
        title: "RTEA Aperture 3 Hexagonal",
        description: "A Discrete Global Grid Reference System based on the Rhombic \
                      Triacontahedron Equal Area projection, with aperture 3 hexagonal zones, \
                      using an indexing scheme based on RTEA9R.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/RTEA3H",
        registered: false,
        crs: None,
        grid: || AnyGrid::Rtea3h(rs4dggs::rtea3h()),
        definition: include_str!("../definitions/RTEA3H.json"),
    },
    Dggrs {
        id: "ISEA7H_Z7",
        title: "ISEA Aperture 7 Hexagonal / Z7 Indexing",
        description: "A Discrete Global Grid Reference System based on the Icosahedral Snyder \
                      Equal Area projection, with aperture 7 hexagonal zones, using the Z7 \
                      indexing scheme.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/ISEA7H_Z7",
        registered: false,
        crs: Some(ISEA_CRS),
        grid: || AnyGrid::Igeo7(rs4dggs::igeo7()),
        definition: include_str!("../definitions/ISEA7H_Z7.json"),
    },
    Dggrs {
        id: "IVEA7H_Z7",
        title: "IVEA Aperture 7 Hexagonal / Z7 Indexing",
        description: "A Discrete Global Grid Reference System based on the Icosahedral \
                      Vertex-oriented Equal Area projection, with aperture 7 hexagonal zones, \
                      using the Z7 indexing scheme.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/IVEA7H_Z7",
        registered: false,
        crs: None,
        grid: || AnyGrid::Ivea7h(rs4dggs::ivea7h()),
        definition: include_str!("../definitions/IVEA7H_Z7.json"),
    },
    Dggrs {
        id: "RTEA7H_Z7",
        title: "RTEA Aperture 7 Hexagonal / Z7 Indexing",
        description: "A Discrete Global Grid Reference System based on the Rhombic \
                      Triacontahedron Equal Area projection, with aperture 7 hexagonal zones, \
                      using the Z7 indexing scheme.",
        uri: "https://www.opengis.net/def/dggrs/OGC/1.0/RTEA7H_Z7",
        registered: false,
        crs: None,
        grid: || AnyGrid::Rtea7h(rs4dggs::rtea7h()),
        definition: include_str!("../definitions/RTEA7H_Z7.json"),
    },
];

/// The longest identifier of the six, `ISEA7H_Z7`, in bytes: a longer text is no identifier.
const LONGEST_ID: usize = 9;

impl Dggrs {
    /// The six grids, in the order ISEA3H, IVEA3H, RTEA3H, ISEA7H_Z7, IVEA7H_Z7, RTEA7H_Z7.
    pub fn all() -> &'static [Dggrs] {
        &ALL
    }

    /// The grid of this identifier, matched exactly (case included), or `None`. `IGEO7` is
    /// accepted as a second name of `ISEA7H_Z7`; `ISEA7H`, `IVEA7H` and `RTEA7H`, which name
    /// grids of another indexing, are not. A text longer than every identifier is refused
    /// before it is compared.
    pub fn from_id(id: &str) -> Option<&'static Dggrs> {
        if id.len() > LONGEST_ID {
            return None;
        }
        let id = if id == "IGEO7" { "ISEA7H_Z7" } else { id };
        ALL.iter().find(|d| d.id == id)
    }

    /// The identifier, as the path of the service carries it: `ISEA7H_Z7` also when the grid
    /// was asked for as `IGEO7`.
    pub fn id(&self) -> &'static str {
        self.id
    }

    /// The title.
    pub fn title(&self) -> &'static str {
        self.title
    }

    /// A description in one sentence.
    pub fn description(&self) -> &'static str {
        self.description
    }

    /// The URI that a data document carries in its member `dggrs`:
    /// `https://www.opengis.net/def/dggrs/OGC/1.0/` followed by the identifier, for all six.
    ///
    /// Only ISEA3H and IVEA3H are registered with the OGC under this URI; for the other four it
    /// is the form the register would give them. Its last segment is the name DGGAL knows the
    /// grid by, which is what DGGAL's own reader of such documents requires.
    pub fn uri(&self) -> &'static str {
        self.uri
    }

    /// The URI under which the OGC has registered the grid: `Some` for ISEA3H and IVEA3H alone,
    /// and then equal to [`Dggrs::uri`]. What a grid's description carries as its `uri`.
    pub fn registered_uri(&self) -> Option<&'static str> {
        self.registered.then_some(self.uri)
    }

    /// The coordinate reference system of the grid's projection, where one is registered: `Some`
    /// for ISEA3H and ISEA7H_Z7 alone.
    pub fn crs(&self) -> Option<&'static str> {
        self.crs
    }

    /// The library's grid.
    pub fn grid(&self) -> AnyGrid {
        (self.grid)()
    }

    /// The definition document of the grid, as JSON text.
    ///
    /// Those of ISEA3H and IVEA3H are the OGC register's files, verbatim. The other four were
    /// written for this crate from the registered definitions of IVEA3H, ISEA7H and IVEA7H and
    /// from DGGAL; each says in its own description that it is not registered with the OGC.
    pub fn definition(&self) -> &'static str {
        self.definition
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_grid_shows_its_identifier_and_not_its_definition() {
        let shown = format!("{:?}", Dggrs::from_id("ISEA3H").unwrap());
        assert_eq!(shown, "Dggrs { id: \"ISEA3H\", .. }");
    }

    #[test]
    fn the_six_are_in_their_order() {
        let ids: Vec<&str> = Dggrs::all().iter().map(Dggrs::id).collect();
        assert_eq!(
            ids,
            [
                "ISEA3H",
                "IVEA3H",
                "RTEA3H",
                "ISEA7H_Z7",
                "IVEA7H_Z7",
                "RTEA7H_Z7"
            ]
        );
    }

    #[test]
    fn each_identifier_names_its_own_entry_and_igeo7_names_isea7h_z7() {
        for dggrs in Dggrs::all() {
            assert!(std::ptr::eq(Dggrs::from_id(dggrs.id()).unwrap(), dggrs));
        }
        assert_eq!(Dggrs::from_id("IGEO7").unwrap().id(), "ISEA7H_Z7");
    }

    #[test]
    fn other_names_and_other_casings_are_refused() {
        let long = "A".repeat(1_000_000);
        for id in [
            "IVEA7H", "RTEA7H", "ISEA7H", "isea3h", "igeo7", "", "ISEA3H ", " ISEA3H", &long,
        ] {
            assert!(Dggrs::from_id(id).is_none(), "{id:.20}");
        }
    }

    #[test]
    fn each_is_the_librarys_grid_of_the_table() {
        let names: Vec<&str> = Dggrs::all().iter().map(|d| d.grid().name()).collect();
        assert_eq!(
            names,
            ["ISEA3H", "IVEA3H", "RTEA3H", "IGEO7", "IVEA7H", "RTEA7H"]
        );
    }

    #[test]
    fn the_uris_end_in_the_identifier_and_two_are_registered() {
        for dggrs in Dggrs::all() {
            assert!(
                dggrs.uri().ends_with(&format!("/{}", dggrs.id())),
                "{}",
                dggrs.uri()
            );
            assert!(
                dggrs
                    .uri()
                    .starts_with("https://www.opengis.net/def/dggrs/OGC/1.0/")
            );
        }
        let registered: Vec<Option<&str>> =
            Dggrs::all().iter().map(Dggrs::registered_uri).collect();
        assert_eq!(
            registered,
            [
                Some("https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H"),
                Some("https://www.opengis.net/def/dggrs/OGC/1.0/IVEA3H"),
                None,
                None,
                None,
                None
            ]
        );
        for dggrs in &Dggrs::all()[..2] {
            assert_eq!(dggrs.registered_uri(), Some(dggrs.uri()));
        }
    }

    #[test]
    fn the_two_isea_grids_carry_a_crs() {
        let crs: Vec<Option<&str>> = Dggrs::all().iter().map(Dggrs::crs).collect();
        let isea = Some("https://www.opengis.net/def/crs/OGC/0/1534");
        assert_eq!(crs, [isea, None, None, isea, None, None]);
    }

    #[test]
    fn no_two_share_an_identifier_a_uri_or_a_grid() {
        let all = Dggrs::all();
        for (i, a) in all.iter().enumerate() {
            for b in &all[i + 1..] {
                assert_ne!(a.id(), b.id());
                assert_ne!(a.uri(), b.uri());
                assert_ne!(a.grid().name(), b.grid().name());
                assert_ne!(a.title(), b.title());
                assert_ne!(a.description(), b.description());
            }
        }
    }

    #[test]
    fn the_registered_definitions_are_the_registers_files() {
        assert_eq!(Dggrs::all()[0].definition().len(), 4880);
        assert_eq!(Dggrs::all()[1].definition().len(), 5037);
    }

    #[test]
    fn each_definition_holds_the_required_members_and_its_own_uri() {
        for dggrs in Dggrs::all() {
            let text = dggrs.definition();
            for member in [
                "\"dggh\"",
                "\"zirs\"",
                "\"textZIRS\"",
                "\"subZoneOrder\"",
                "\"spatialDimensions\"",
                "\"temporalDimensions\"",
            ] {
                assert!(text.contains(member), "{}: {member}", dggrs.id());
            }
            assert!(
                text.contains(&format!("\"uri\": \"{}\"", dggrs.uri())),
                "{}",
                dggrs.id()
            );
            assert!(
                text.contains(&format!(
                    "\"refinementRatio\": {}",
                    dggrs.grid().refinement_ratio()
                )),
                "{}",
                dggrs.id()
            );
            assert_eq!(
                text.contains("not registered with the OGC"),
                dggrs.registered_uri().is_none(),
                "{}",
                dggrs.id()
            );
        }
    }
}
