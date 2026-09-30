//! CSV rows: a header, then fields joined by commas. The fields are numbers and zone identifiers,
//! which hold no comma, quotation mark or line break, so nothing is quoted.

use std::io::{self, Write};

pub struct Csv<'a> {
    out: &'a mut dyn Write,
    header: &'static str,
    started: bool,
}

impl<'a> Csv<'a> {
    pub fn new(out: &'a mut dyn Write, header: &'static str) -> Self {
        Csv {
            out,
            header,
            started: false,
        }
    }

    /// Writes the header, once. It goes out with the first row, or at the end, so that an input
    /// refused outright writes nothing.
    pub fn begin(&mut self) -> io::Result<()> {
        if !self.started {
            self.started = true;
            writeln!(self.out, "{}", self.header)?;
        }
        Ok(())
    }

    pub fn row(&mut self, fields: &[&str]) -> io::Result<()> {
        self.begin()?;
        writeln!(self.out, "{}", fields.join(","))
    }
}
