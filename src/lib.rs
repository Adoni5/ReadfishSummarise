#![deny(missing_docs)]
#![warn(clippy::missing_docs_in_private_items)]
#![allow(dead_code)]
//! # readfish_summarise
//!
//! `readfish_summarise` is a collection of utilities to provide a standardised way of summarising
//! readfish runs that have been run. Currently the accepted analysable inputs are PAF records from running the readfish
//! summarise entrypoint.
//!
//! The crate is split into modules handling separate functionalities.
//!
use std::{cell::RefCell, collections::HashMap, error::Error, fmt, ops::Deref};

use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use prettytable::{color, row, Attr, Cell, Row, Table};
use pyo3::{prelude::*, types::PyIterator};
use readfish_tools::nanopore::{format_bases, running_mean};
use readfish_tools::paf::PafRecord;

/// Dynamic result type for holding either a generic value or an error
pub type DynResult<T> = Result<T, Box<dyn Error + 'static>>;

/// Colour of the on target columns in table
const FG_ON: Attr = Attr::ForegroundColor(color::BRIGHT_CYAN);
/// COlour of the off target columns in table
const FG_OFF: Attr = Attr::ForegroundColor(color::BRIGHT_MAGENTA);
/// Colour of other columns
const FG_OTHER: Attr = Attr::ForegroundColor(color::BRIGHT_WHITE);
/// Represents the mean read lengths for on-target, off-target, and total reads.
#[derive(Debug)]
pub struct MeanReadLengths {
    /// The mean read length of on-target reads.
    pub on_target: isize,
    /// Number of on target reads analysed
    on_target_count: isize,
    /// The mean read length of off-target reads.
    pub off_target: isize,
    /// Number of off target reads analysed
    off_target_count: isize,
    /// The mean read length of all reads (on-target + off-target).
    pub total: isize,
    /// Number of reads analysed
    total_count: isize,
}

impl MeanReadLengths {
    /// Creates a new `MeanReadLengths` instance with all fields initialized to 0.
    pub fn new() -> Self {
        MeanReadLengths {
            on_target: 0,
            on_target_count: 0,
            off_target: 0,
            off_target_count: 0,
            total: 0,
            total_count: 0,
        }
    }

    /// Updates the mean read lengths for on-target, off-target, and total reads based on the provided
    /// PAF record and whether the read is on-target or off-target.
    ///
    /// # Arguments
    ///
    /// * `paf` - A reference to the [`PafRecord`] representing the alignment record for a read.
    /// * `on_target` - A boolean indicating whether the read is on-target (true) or off-target (false).
    ///
    /// # Example
    ///
    /// ```
    /// use readfish_tools::{MeanReadLengths, paf::PafRecord};
    /// let mut mean_lengths = MeanReadLengths::new();
    /// let paf_record = PafRecord::new("read123 200 0 200 + contig123 300 0 300 200 200 50 ch=1".split(" ").collect()).unwrap();
    /// mean_lengths.update_lengths(&paf_record, true);
    /// ```
    pub fn update_lengths(&mut self, paf: &PafRecord, on_target: bool) {
        if on_target {
            running_mean(
                &mut self.on_target,
                &mut self.on_target_count,
                &mut (paf.query_length as isize),
            );
        } else {
            running_mean(
                &mut self.off_target,
                &mut self.off_target_count,
                &mut (paf.query_length as isize),
            );
        }
        running_mean(
            &mut self.total,
            &mut self.total_count,
            &mut (paf.query_length as isize),
        );
    }
}

impl Default for MeanReadLengths {
    fn default() -> Self {
        Self::new()
    }
}
/// Divide too floats, returning none if the divisor is zero
fn safe_divide(dividend: f64, divisor: f64) -> Option<f64> {
    if divisor == 0.0 {
        None // Division by zero, return None or handle it as needed
    } else {
        Some(dividend / divisor)
    }
}

/// Represents a summary of a contig or sequence from a sequencing experiment.
/// It includes various metrics related to the contig's characteristics and read mapping.
#[derive(Debug)]
pub struct ContigSummary {
    /// The name or identifier of the contig.
    pub name: String,
    /// The length of the contig in base pairs.
    pub length: usize,
    /// Number of reads total that mapped
    pub total_reads: usize,
    /// The mean read length of the mapped reads associated with this contig.
    pub mean_read_lengths: MeanReadLengths,
    /// The count of reads that are mapped on the target region (on-target reads).
    pub on_target_alignment_count: usize,
    /// The count of reads that are mapped off the target region (off-target reads).
    pub off_target_alignment_count: usize,
    /// The total yield (base pairs) of on-target reads for this contig.
    pub off_target_yield: usize,
    /// The total yield (base pairs) of off-target reads for this contig.
    pub on_target_yield: usize,
    /// number of targets on contig
    pub number_of_targets: usize,
    /// Number of bases in target regions, used to calculate rough on target coverage
    pub number_target_bases: usize,
}
impl ContigSummary {
    /// Create a new `ContigSummary` instance with default values for all fields except `name` and `length`.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the contig.
    /// * `length` - The length of the contig.
    pub fn new(name: String, length: usize) -> Self {
        ContigSummary {
            name,
            length,
            mean_read_lengths: MeanReadLengths::new(),
            total_reads: 0,
            on_target_alignment_count: 0,
            off_target_alignment_count: 0,
            on_target_yield: 0,
            off_target_yield: 0,
            number_of_targets: 0,
            number_target_bases: 0,
        }
    }
    /// Total number of alignments
    pub fn total_alignments(&self) -> usize {
        self.off_target_alignment_count + self.on_target_alignment_count
    }

    /// Percent of on target alignments against total alignments
    pub fn on_target_alignment_percent(&self) -> f64 {
        self.on_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }

    /// Percent of off target alignments against total alignments
    pub fn off_target_alignment_percent(&self) -> f64 {
        self.off_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }
    /// Percent of on target yield against total yield
    pub fn on_target_yield_percent(&self) -> f64 {
        self.on_target_yield as f64 / self.total_yield() as f64 * 100.0
    }

    /// Percent of off target yield against total yield
    pub fn off_target_yield_percent(&self) -> f64 {
        self.off_target_yield as f64 / self.total_yield() as f64 * 100.0
    }
    /// Total yield
    pub fn total_yield(&self) -> usize {
        self.on_target_yield + self.off_target_yield
    }
    /// Yield ratio (on target:off target)
    pub fn yield_ratio(&self) -> String {
        format!(
            "{}:{:.2}",
            safe_divide(self.on_target_yield as f64, self.on_target_yield as f64).unwrap_or(0.0),
            safe_divide(self.off_target_yield as f64, self.on_target_yield as f64).unwrap_or(0.0)
        )
        .to_string()
    }

    /// Estimated on target coverage
    pub fn estimated_target_coverage(&self) -> String {
        let yieldy = self
            .on_target_yield
            .checked_div(self.number_target_bases)
            .unwrap_or(0);
        format_bases(yieldy)
    }
    /// Mean read length of all reads on the contig.
    pub fn mean_read_length(&self) -> usize {
        self.mean_read_lengths.total as usize
    }
    /// On target mean read length of all reads on the contig.
    pub fn on_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths.on_target as usize
    }
    /// Off target mean read length of all reads on the contig.
    pub fn off_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths.off_target as usize
    }

    /// Add a readfish target to the Contig summary
    pub fn add_target(&mut self, start: usize, end: usize) -> DynResult<()> {
        self.number_of_targets += 1;
        self.number_target_bases += end - start;
        Ok(())
    }

    /// Add a given number of total reads
    pub fn add_total_reads(&mut self, total_reads: usize) {
        self.total_reads += total_reads;
    }
}
#[derive(Debug)]
/// Represents a summary of sequencing data, including various metrics related to the output of the experiment.
pub struct ConditionSummary {
    /// The name or identifier of the sequencing data.
    pub name: String,
    /// The total number of reads in the sequencing data.
    pub total_reads: usize,
    /// Mean read lengths
    pub mean_read_lengths: MeanReadLengths,
    /// The count of reads that are mapped off the target regions (off-target reads).
    pub off_target_alignment_count: usize,
    /// The count of reads that are mapped to the target regions (on-target reads).
    pub on_target_alignment_count: usize,
    /// The total yield (base pairs) of off-target reads in the sequencing data.
    pub off_target_yield: usize,
    /// The total yield (base pairs) of on-target reads in the sequencing data.
    pub on_target_yield: usize,
    /// The mean read quality of off-target reads.
    pub off_target_mean_read_quality: f64,
    /// The mean read quality of on-target reads.
    pub on_target_mean_read_quality: f64,
    /// The N50 metric for the entire dataset, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length.
    pub n50: usize,
    /// The N50 metric for on-target reads, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length for on-target reads.
    pub on_target_n50: usize,
    /// The N50 metric for off-target reads, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length for off-target reads.
    pub off_target_n50: usize,
    /// Number of target regions
    pub number_of_targets: usize,
    /// Number of bases covered by targets
    pub number_target_bases: usize,
    /// A vector of `ContigSummary` representing summaries of individual contigs or sequences
    /// in the sequencing data.
    pub contigs: HashMap<String, ContigSummary>,
}

impl fmt::Display for ConditionSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Condition Name: {}", self.name)?;
        writeln!(f, "Total Reads: {}", self.total_reads)?;
        writeln!(
            f,
            "Off-Target Read Count: {}",
            self.off_target_alignment_count
        )?;
        writeln!(
            f,
            "On-Target Read Count: {}",
            self.on_target_alignment_count
        )?;
        writeln!(
            f,
            "Off-Target Percent: {:.2}%",
            self.off_target_yield_percent()
        )?;
        writeln!(f, "Off-Target Yield: {}", self.off_target_yield)?;
        writeln!(f, "On-Target Yield: {}", self.on_target_yield)?;
        writeln!(
            f,
            "Off-Target Mean Read Length: {}",
            self.off_target_mean_read_length()
        )?;
        writeln!(
            f,
            "On-Target Mean Read Length: {}",
            self.on_target_mean_read_length()
        )?;
        // writeln!(
        //     f,
        //     "Off-Target Mean Read Quality: {:.2}",
        //     self.off_target_mean_read_quality
        // )?;
        // writeln!(
        //     f,
        //     "On-Target Mean Read Quality: {:.2}",
        //     self.on_target_mean_read_quality
        // )?;
        // writeln!(f, "N50: {}", self.n50)?;
        // writeln!(f, "On-Target N50: {}", self.on_target_n50)?;
        // writeln!(f, "Off-Target N50: {}", self.off_target_n50)?;

        writeln!(f, "Contigs:")?;
        for (contig_name, contig_summary) in &self.contigs {
            writeln!(f, "  Contig Name: {}", contig_name)?;
            writeln!(f, "  Length: {}", contig_summary.length)?;
            // Print other fields from ContigSummary here
            // For example:
            // writeln!(f, "  Contig Mean Read Length: {}", contig_summary.mean_read_length)?;
        }
        Ok(())
    }
}

impl ConditionSummary {
    /// Add a target to the COnditions summary and the correct contig summary
    pub fn add_target(
        &mut self,
        contig: String,
        contig_len: usize,
        start: usize,
        end: usize,
    ) -> DynResult<()> {
        self.number_of_targets += 1;
        self.number_target_bases += end - start;
        let contig = self.get_or_add_contig(&contig, contig_len);
        contig.add_target(start, end)?;
        Ok(())
    }

    /// Update the `ConditionSummary` with information from the provided `PafRecord`.
    ///
    /// This method updates the fields of the `ConditionSummary` based on the information
    /// from the given `PafRecord`. It increments the appropriate read counts (on-target
    /// or off-target), calculates the mean read lengths and read qualities, updates the
    /// total reads count, and calculates the off-target percentage.
    ///
    /// # Arguments
    ///
    /// * `paf` - The [`PafRecord`] containing the information about the alignment.
    /// * `on_target` - A boolean flag indicating whether the alignment is on-target or off-target.
    ///
    /// # Returns
    ///
    /// This function returns a [`DynResult`] (a dynamic result that can contain any error).
    /// If the operation is successful, the `DynResult` will hold an `Ok(())`. Otherwise, it
    /// will hold an `Err` containing a helpful error message.
    pub fn update(&mut self, paf: PafRecord, on_target: bool) -> DynResult<()> {
        // update the condition struct
        self.mean_read_lengths.update_lengths(&paf, on_target);
        if on_target {
            self.on_target_alignment_count += 1;
            self.on_target_yield += paf.query_length;
            // self.on_target_mean_read_quality += paf.tlen as f64;
        } else {
            self.off_target_alignment_count += 1;
            self.off_target_yield += paf.query_length;
            // self.off_target_mean_read_quality += paf.tlen as f64;
        }
        let contig = self.get_or_add_contig(&paf.target_name, paf.target_length);
        contig.mean_read_lengths.update_lengths(&paf, on_target);
        if on_target {
            contig.on_target_alignment_count += 1;
            contig.on_target_yield += paf.query_length;
            // self.on_target_mean_read_quality += paf.tlen as f64;
        } else {
            contig.off_target_alignment_count += 1;
            contig.off_target_yield += paf.query_length;
            // self.off_target_mean_read_quality += paf.tlen as f64;
        }
        // contig.mean_read_quality = paf.tlen;
        // contig.n50 = paf.tlen;
        // contig.on_target_read_count = paf.tlen;
        // contig.off_target_read_count = paf.tlen;

        Ok(())
    }
    /// Create a new `Summary` instance with default values for all fields except `name`.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the summary.
    pub fn new(name: String) -> Self {
        ConditionSummary {
            name,
            total_reads: 0,
            off_target_alignment_count: 0,
            on_target_alignment_count: 0,
            off_target_yield: 0,
            on_target_yield: 0,
            mean_read_lengths: MeanReadLengths::new(),
            off_target_mean_read_quality: 0.0,
            on_target_mean_read_quality: 0.0,
            n50: 0,
            on_target_n50: 0,
            off_target_n50: 0,
            number_of_targets: 0,
            number_target_bases: 0,
            contigs: HashMap::new(),
        }
    }
    /// The percentage of off-target reads in the sequencing data.
    pub fn off_target_yield_percent(&self) -> f64 {
        self.off_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }

    /// The percentage of off-target reads in the sequencing data.
    pub fn on_target_yield_percent(&self) -> f64 {
        self.on_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }

    /// Get the name or identifier of the sequencing data.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the name or identifier of the sequencing data.
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Get the total number of reads in the sequencing data.
    pub fn total_reads(&self) -> usize {
        self.total_reads
    }

    /// Total number of alignments (on and off target combined)
    pub fn total_alignments(&self) -> usize {
        self.off_target_alignment_count + self.on_target_alignment_count
    }

    /// Percent of alignments that were off target
    pub fn off_target_alignment_percent(&self) -> f64 {
        self.off_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }

    /// Percent of alignments that were on target
    pub fn on_target_alignment_percent(&self) -> f64 {
        self.on_target_alignment_count as f64 / self.total_alignments() as f64 * 100.0
    }

    /// Set the total number of reads in the sequencing data.
    pub fn add_total_reads(&mut self, total_reads: usize) {
        self.total_reads += total_reads;
    }

    /// Get the count of reads that are mapped off the target regions (off-target reads).
    pub fn off_target_read_count(&self) -> usize {
        self.off_target_alignment_count
    }

    /// Set the count of reads that are mapped off the target regions (off-target reads).
    pub fn set_off_target_read_count(&mut self, off_target_read_count: usize) {
        self.off_target_alignment_count = off_target_read_count;
    }

    /// Get the count of reads that are mapped to the target regions (on-target reads).
    pub fn on_target_read_count(&self) -> usize {
        self.on_target_alignment_count
    }

    /// Set the count of reads that are mapped to the target regions (on-target reads).
    pub fn set_on_target_read_count(&mut self, on_target_read_count: usize) {
        self.on_target_alignment_count = on_target_read_count;
    }

    /// Get the total yield (base pairs) of off-target reads in the sequencing data.
    pub fn off_target_yield(&self) -> usize {
        self.off_target_yield
    }

    /// Set the total yield (base pairs) of off-target reads in the sequencing data.
    pub fn set_off_target_yield(&mut self, off_target_yield: usize) {
        self.off_target_yield = off_target_yield;
    }

    /// Get the total yield (base pairs) of on-target reads in the sequencing data.
    pub fn on_target_yield(&self) -> usize {
        self.on_target_yield
    }

    /// Set the total yield (base pairs) of on-target reads in the sequencing data.
    pub fn set_on_target_yield(&mut self, on_target_yield: usize) {
        self.on_target_yield = on_target_yield;
    }
    /// Get the mean read length of all reads
    pub fn mean_read_length(&self) -> usize {
        self.mean_read_lengths.total as usize
    }

    /// Get the mean read length of off-target reads.
    pub fn off_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths.off_target as usize
    }

    /// Get the mean read length of on-target reads.
    pub fn on_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths.on_target as usize
    }

    /// Get the mean read quality of off-target reads.
    pub fn off_target_mean_read_quality(&self) -> f64 {
        self.off_target_mean_read_quality
    }

    /// Set the mean read quality of off-target reads.
    pub fn set_off_target_mean_read_quality(&mut self, off_target_mean_read_quality: f64) {
        self.off_target_mean_read_quality = off_target_mean_read_quality;
    }

    /// Get the mean read quality of on-target reads.
    pub fn on_target_mean_read_quality(&self) -> f64 {
        self.on_target_mean_read_quality
    }

    /// Set the mean read quality of on-target reads.
    pub fn set_on_target_mean_read_quality(&mut self, on_target_mean_read_quality: f64) {
        self.on_target_mean_read_quality = on_target_mean_read_quality;
    }

    /// Get the N50 metric for the entire dataset.
    pub fn n50(&self) -> usize {
        self.n50
    }

    /// Set the N50 metric for the entire dataset.
    pub fn set_n50(&mut self, n50: usize) {
        self.n50 = n50;
    }

    /// Get the N50 metric for on-target reads.
    pub fn on_target_n50(&self) -> usize {
        self.on_target_n50
    }

    /// Set the N50 metric for on-target reads.
    pub fn set_on_target_n50(&mut self, on_target_n50: usize) {
        self.on_target_n50 = on_target_n50;
    }

    /// Get the N50 metric for off-target reads.
    pub fn off_target_n50(&self) -> usize {
        self.off_target_n50
    }

    /// Set the N50 metric for off-target reads.
    pub fn set_off_target_n50(&mut self, off_target_n50: usize) {
        self.off_target_n50 = off_target_n50;
    }

    /// Get a reference to the vector of `ContigSummary`.
    pub fn contigs(&self) -> &HashMap<String, ContigSummary> {
        &self.contigs
    }

    /// Get a mutable reference to the vector of `ContigSummary`.
    pub fn contigs_mut(&mut self) -> &mut HashMap<String, ContigSummary> {
        &mut self.contigs
    }
    /// Get the on_target/off target ratio
    pub fn yield_ratio(&self) -> String {
        format!(
            "{}:{}",
            safe_divide(self.on_target_yield as f64, self.on_target_yield as f64).unwrap_or(0.0),
            safe_divide(self.off_target_yield as f64, self.on_target_yield as f64).unwrap_or(0.0)
        )
        .to_string()
    }

    /// Estimate on target coverage
    pub fn estimated_target_coverage(&self) -> String {
        let yieldy = self
            .on_target_yield()
            .checked_div(self.number_target_bases)
            .unwrap_or(0);
        format_bases(yieldy)
    }

    /// Get the ContigSummary associated with the given contig name or
    ///  add a new ContigSummary with the specified name and length if it doesn't exist.
    ///
    /// # Arguments
    ///
    /// * `contig` - The name of the contig.
    /// * `length` - The length of the contig.
    ///
    /// # Returns
    ///
    /// A reference to the ContigSummary associated with the contig name.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::collections::HashMap;
    /// use your_crate::ContigSummary;
    ///
    /// let mut contig_map: HashMap<String, ContigSummary> = HashMap::new();
    ///
    /// // Get an existing contig or add a new one if it doesn't exist
    /// let contig_name = "chr1".to_string();
    /// let contig_length = 1000;
    /// let contig_summary = contig_map.get_or_add_contig(contig_name.clone(), contig_length);
    ///
    /// // Now you can modify the ContigSummary fields or access its properties
    /// println!("Contig name: {}", contig_summary.name);
    /// println!("Contig length: {}", contig_summary.length);
    /// ```
    pub fn get_or_add_contig(&mut self, contig: &str, length: usize) -> &mut ContigSummary {
        self.contigs
            .entry(contig.to_string())
            .or_insert(ContigSummary::new(contig.to_string(), length))
    }

    /// get the total yield
    pub fn total_yield(&self) -> usize {
        self.on_target_yield + self.off_target_yield
    }
}

/// A struct representing a summary of conditions.
///
/// The `Summary` struct contains a hashmap where each key represents the name of a condition, and the corresponding value is a `ConditionSummary` struct
/// containing the summary information for that condition.
///
/// # Fields
///
/// * `conditions` - A hashmap containing the summary information for each condition. The key is a string representing the name of the condition,
/// and the value is a `ConditionSummary` struct containing the summary information for that condition.
///
/// # Examples
///
/// ```rust, ignore
/// use std::collections::HashMap;
/// use readfish_tools::{Summary, ConditionSummary};
///
/// // Create a new Summary
/// let mut summary = Summary {
///     conditions: HashMap::new(),
/// };
///
/// // Add some condition summaries
/// summary.conditions.insert(
///     "ConditionA".to_string(),
///     ConditionSummary {
///         // ... fill in the details for ConditionA ...
///     }
/// );
///
/// summary.conditions.insert(
///     "ConditionB".to_string(),
///     ConditionSummary {
///         // ... fill in the details for ConditionB ...
///     }
/// );
///
/// // Access a specific condition summary
/// if let Some(condition_summary) = summary.conditions.get("ConditionA") {
///     println!("Summary for ConditionA: {:?}", condition_summary);
/// }
/// ```
#[derive(Debug)]
pub struct Summary {
    /// Conditions summary for a given region or barcode.
    pub conditions: HashMap<String, ConditionSummary>,
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Todo rewrite to use Macro!
        let mut condition_table = Table::new();
        condition_table.add_row(Row::new(vec![
            Cell::new("Condition")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Reads")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Alignments")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Yield")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(4),
            Cell::new("Mean read\n lengths")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Number of\n targets")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Estimated\n coverage")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
        ]));
        condition_table.add_row(row!["", "", FBb->"On-Target", FMb->"Off-Target", FWb->"Total",  FBb->"On-Target", FMb->"Off-Target", FWb->"Total", FWb->"Ratio", FBb->"On-target", FMb->"Off-target", FWb->"Combined"]);
        for (condition_name, condition_summary) in &self.conditions {
            condition_table.add_row(Row::new(vec![
                Cell::new(condition_name).with_style(Attr::ForegroundColor(color::BRIGHT_YELLOW)),
                // total reads
                Cell::new(
                    &condition_summary
                        .total_reads
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OTHER),
                // on target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .on_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    condition_summary.on_target_alignment_percent()
                ))
                .with_style(FG_ON),
                // off target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .off_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    condition_summary.off_target_alignment_percent()
                ))
                .with_style(FG_OFF),
                // total alignments
                Cell::new(
                    &condition_summary
                        .total_alignments()
                        .to_formatted_string(&Locale::en)
                        .to_string(),
                )
                .with_style(FG_OTHER),
                // on target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    &format_bases(condition_summary.on_target_yield),
                    condition_summary.on_target_yield_percent()
                ))
                .with_style(FG_ON),
                // off target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    &format_bases(condition_summary.off_target_yield),
                    condition_summary.off_target_yield_percent()
                ))
                .with_style(FG_OFF),
                // total yield
                Cell::new(&format_bases(condition_summary.total_yield())).with_style(FG_OTHER),
                Cell::new(&condition_summary.yield_ratio()).with_style(FG_OTHER),
                // on target mean read length
                Cell::new(&format_bases(
                    condition_summary.on_target_mean_read_length(),
                ))
                .with_style(FG_ON),
                // off target mean read length
                Cell::new(&format_bases(
                    condition_summary.off_target_mean_read_length(),
                ))
                .with_style(FG_OFF),
                // mean read length
                Cell::new(&format_bases(condition_summary.mean_read_length())).with_style(FG_OTHER),
                // number of targets
                Cell::new(
                    &condition_summary
                        .number_of_targets
                        .to_formatted_string(&Locale::en)
                        .to_string(),
                )
                .with_style(FG_OTHER),
                // Estimated target coverage
                Cell::new(&condition_summary.estimated_target_coverage()).with_style(FG_OTHER),
            ]));

            // writeln!(
            //     f,
            //     "  Off-Target Mean Read Quality: {:.2}",
            //     condition_summary.off_target_mean_read_quality
            // )?;
            // writeln!(
            //     f,
            //     "  On-Target Mean Read Quality: {:.2}",
            //     condition_summary.on_target_mean_read_quality
            // )?;
            // writeln!(f, "  N50: {}", condition_summary.n50)?;
            // writeln!(f, "  On-Target N50: {}", condition_summary.on_target_n50)?;
            // writeln!(f, "  Off-Target N50: {}", condition_summary.off_target_n50)?;
        }
        condition_table.add_row(row!["", "", FBb->"On-Target", FMb->"Off-Target", FWb->"Total",  FBb->"On-Target", FMb->"Off-Target", FWb->"Total", FWb->"Ratio", FBb->"On-target", FMb->"Off-target", FWb->"Combined"]);
        condition_table.add_row(Row::new(vec![
            Cell::new("Condition")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Reads")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Alignments")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Yield")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(4),
            Cell::new("Mean read\n lengths")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Number of\n targets")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Estimated\n coverage")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
        ]));
        condition_table.printstd();
        writeln!(f, "Contigs:")?;

        for condition_summary in self.conditions.values() {
            let mut contig_table = Table::new();
            contig_table.add_row(Row::new(vec![
                Cell::new("Condition Name")
                    .with_style(Attr::Bold)
                    .with_style(FG_OTHER),
                Cell::new(&condition_summary.name)
                    .with_style(Attr::Standout(true))
                    .with_style(FG_OTHER)
                    .with_hspan(14),
            ]));
            // Create a custom format with left-leading spaces
            contig_table.get_format();
            contig_table.add_row(Row::new(vec![
                Cell::new("Contig")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
                Cell::new("Contig Length")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
                Cell::new("Reads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
                Cell::new("Alignments")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                    .with_hspan(3),
                Cell::new("Yield")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                    .with_hspan(4),
                Cell::new("Mean \nRead Length")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                    .with_hspan(3),
                Cell::new("Number of\n targets")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
                Cell::new("Estimated\n coverage")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            ]));
            contig_table.add_row(row!["", "", "", FBb->"On-Target", FMb->"Off-Target", b->"Total",  FBb->"On-Target", FMb->"Off-Target", b->"Total", b->"Ratio", FBb->"On-target", FMb->"Off-target", b->"Combined"]);

            for (contig_name, contig_summary) in condition_summary
                .contigs
                .iter()
                .sorted_by(|(key1, _), (key2, _)| natord::compare(key1, key2))
            {
                contig_table.add_row(Row::new(vec![
                    // contig name
                    Cell::new(contig_name)
                        .with_style(Attr::Blink)
                        .with_style(FG_OTHER),
                    // contig length
                    Cell::new(&contig_summary.length.to_formatted_string(&Locale::en))
                        .with_style(FG_OTHER),
                    // Number of reads
                    Cell::new(&contig_summary.total_reads.to_formatted_string(&Locale::en))
                        .with_style(FG_OTHER),
                    // on target alignment
                    Cell::new(&format!(
                        "{} ({:.2}%)",
                        contig_summary
                            .on_target_alignment_count
                            .to_formatted_string(&Locale::en),
                        contig_summary.on_target_alignment_percent()
                    ))
                    .with_style(FG_ON),
                    // off target alignment
                    Cell::new(&format!(
                        "{} ({:.2}%)",
                        contig_summary
                            .off_target_alignment_count
                            .to_formatted_string(&Locale::en),
                        contig_summary.off_target_alignment_percent()
                    ))
                    .with_style(FG_OFF),
                    // total alignments
                    Cell::new(
                        &contig_summary
                            .total_alignments()
                            .to_formatted_string(&Locale::en),
                    )
                    .with_style(FG_OTHER),
                    // on target yield
                    Cell::new(&format!(
                        "{} ({:.2}%)",
                        contig_summary
                            .on_target_yield
                            .to_formatted_string(&Locale::en),
                        contig_summary.on_target_yield_percent()
                    ))
                    .with_style(FG_ON),
                    //off target yield
                    Cell::new(&format!(
                        "{} ({:.2}%)",
                        contig_summary
                            .off_target_yield
                            .to_formatted_string(&Locale::en),
                        contig_summary.off_target_yield_percent()
                    ))
                    .with_style(FG_OFF),
                    // total yield
                    Cell::new(&format_bases(contig_summary.total_yield())).with_style(FG_OTHER),
                    // yield ratio
                    Cell::new(&contig_summary.yield_ratio()).with_style(FG_OTHER),
                    // on target mean read length
                    Cell::new(&format_bases(contig_summary.on_target_mean_read_length()))
                        .with_style(FG_ON),
                    // off target mean read length
                    Cell::new(&format_bases(contig_summary.off_target_mean_read_length()))
                        .with_style(FG_OFF),
                    // mean read length
                    Cell::new(&format_bases(contig_summary.mean_read_length()))
                        .with_style(FG_OTHER),
                    // number of targets
                    Cell::new(
                        &contig_summary
                            .number_of_targets
                            .to_formatted_string(&Locale::en)
                            .to_string(),
                    )
                    .with_style(FG_OTHER),
                    // estimated target coverage
                    Cell::new(&contig_summary.estimated_target_coverage()).with_style(FG_OTHER),
                ]));
                // Print other fields from ContigSummary here
                // For example:
                // writeln!(f, "    Contig Mean Read Length: {}", contig_summary.mean_read_length)?;
            }
            contig_table.printstd();
        }
        Ok(())
    }
}

impl Summary {
    /// Create a new `Summary` instance with default values for all fields.
    fn new() -> Self {
        Summary {
            conditions: HashMap::new(),
        }
    }

    /// Get the summary for the specified condition. If the condition does not exist in the
    /// `Summary`, it will be created with default values.
    ///
    /// # Arguments
    ///
    /// * `condition_name`: A `String` representing the name or identifier of the condition.
    ///
    /// # Returns
    ///
    /// A reference to the `ConditionSummary` for the specified condition.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::collections::HashMap;
    ///
    /// let mut summary = Summary::new();
    ///
    /// // Get or add the condition with the name "Condition A"
    /// let condition_a = summary.conditions("Condition A".to_string());
    ///
    /// // Modify the fields of the condition summary
    /// condition_a.set_total_reads(10000);
    /// condition_a.set_on_target_read_count(8000);
    /// condition_a.set_off_target_read_count(2000);
    /// // ...
    ///
    /// // Get or add another condition
    /// let condition_b = summary.conditions("Condition B".to_string());
    /// // ...
    /// ```
    pub fn conditions<T: Deref<Target = str>>(
        &mut self,
        condition_name: T,
    ) -> &mut ConditionSummary {
        self.conditions
            .entry(condition_name.to_string())
            .or_insert(ConditionSummary::new(condition_name.to_string()))
    }
}

// PYTHON PyO3 STuff below ////////////////////////

#[pyclass]
#[derive(Debug, Clone)]
/// Stores metadata about a read's mapping and condition.
pub struct MetaData {
    /// The name of the condition from ReadFish analysis.
    pub condition_name: String,
    /// Indicates whether the read mapped to an on-target region of the genome.
    pub on_target: bool,
    /// The Pafline to be analysed
    pub paf_line: String,
}

#[pymethods]
impl MetaData {
    /// Creates a new instance of MetaData with the specified values.
    ///
    /// # Arguments
    ///
    /// * `condition_name` - The name of the condition from ReadFish analysis.
    /// * `on_target` - Indicates whether the read mapped to an on-target region of the genome.
    /// * `paf_line` - Correctly formatted PAF line.
    ///
    /// # Returns
    ///
    /// A new MetaData instance.
    #[new]
    #[pyo3(signature = (condition_name, on_target, paf_line))]
    fn py_new(condition_name: String, on_target: bool, paf_line: String) -> PyResult<Self> {
        Ok(MetaData {
            condition_name,
            on_target,
            paf_line,
        })
    }
}

#[pyclass]
/// Organise the data and methods for analysing a readfish PAF file.
pub struct ReadfishSummary {
    /// Stores the aggregated summary numbers for the readfish run
    summary: RefCell<Summary>,
}

impl Default for ReadfishSummary {
    fn default() -> Self {
        ReadfishSummary::new()
    }
}

impl ReadfishSummary {
    /// Creates a new instance of `ReadfishSummary` with default values.
    ///
    /// This function initializes a new `ReadfishSummary` struct with default values
    /// for all fields. The `summary` field will be initialized with an empty `Summary`
    /// instance. The `_conf`, `_sequencing_summary`, and `_paf_file` fields will be set
    /// to `None`, indicating that they have not been initialized with specific values yet.
    ///
    /// # Returns
    ///
    /// A new `ReadfishSummary` instance with default values.
    ///
    /// # Examples
    ///
    /// ```
    /// use readfish_tools::ReadfishSummary;
    ///
    /// let summary = ReadfishSummary::new();
    /// assert_eq!(summary.has_conf(), false); // _conf field is not set yet
    /// assert_eq!(summary.has_sequencing_summary(), false); // _sequencing_summary field is not set yet
    /// assert_eq!(summary.has_paf_file(), false); // _paf_file field is not set yet
    /// ```
    pub fn new() -> Self {
        ReadfishSummary {
            summary: RefCell::new(Summary::new()),
        }
    }
    // /// Update a condition on the summary
    // pub fn update_condition(
    //     &mut self,
    //     condition_name: &str,
    //     paf_record: PafRecord,
    //     on_target: bool,
    // ) {
    //     let condition_summary = self.summary.borrow_mut().conditions(condition_name);
    //     condition_summary.update(paf_record, on_target).unwrap();
    // }
}

/// Implements methods for interacting with a ReadfishSummary instance from Python.
#[pymethods]
impl ReadfishSummary {
    /// Creates a new instance of ReadfishSummary with default values.
    /// Returns:
    ///     A new ReadfishSummary instance.
    #[new]
    #[pyo3(signature = ())]
    fn py_new() -> PyResult<Self> {
        Ok(ReadfishSummary::default())
    }

    /// Parses PAF lines from a Python iterator and updates the ReadfishSummary accordingly.
    ///
    /// This method takes a Python iterator that provides PAF lines as strings. It iterates over the lines,
    /// parses each line to extract relevant information, and updates the ReadfishSummary based on the
    /// extracted data. The extracted metadata is used to make decisions and update the internal state of
    /// the ReadfishSummary. Finishes by printing a summary of the parsed PAF files to stdout.
    ///
    /// # Arguments
    ///
    /// * `iter`: A Python iterator that provides PAF lines as strings.
    ///
    /// # Returns
    ///
    /// A `PyResult` indicating success or an error encountered during parsing.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # use pyo3::types::PyIterator;
    /// # use readfish_tools::{ReadfishSummary, paf::Metadata, paf::_parse_paf_line};
    ///
    /// // Assuming we have valid inputs
    /// let mut readfish_summary = ReadfishSummary::default();
    /// let paf_lines: Vec<String> = vec![
    ///     "read123 200 0 200 + contig123 300 0 300 200 200 50 ch=1".to_string(),
    ///     "read456 150 0 150 - contig456 200 0 200 150 150 45 ch=2 ba=sampleB".to_string(),
    ///     // Add more PAF lines as needed
    /// ];
    /// let py_iter = PyIterator::new(paf_lines.into_iter());
    ///
    /// let result = readfish_summary.parse_paf_from_iter(&py_iter);
    ///
    /// assert!(result.is_ok());
    /// ```
    fn parse_paf_from_iter(&mut self, iter: &PyIterator) -> PyResult<()> {
        for meta_data in iter {
            let meta_data = meta_data?;
            let meta_data: MetaData = meta_data.extract()?;
            self.update_summary(meta_data)?;
        }
        Ok(())
    }
    /// Updates the summary with information from a given `MetaData`.
    ///
    /// This function updates the summary information based on the provided `MetaData` object,
    /// which contains metadata about a read's mapping and condition. The function extracts
    /// relevant information from the `MetaData` and the associated PAF record, and updates
    /// the summary with the data.
    ///
    /// # Arguments
    ///
    /// * `meta_data` - A `MetaData` object containing information about the read's mapping
    ///   and condition.
    ///
    /// # Returns
    ///
    /// A `PyResult` indicating the success or failure of the update operation.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use readfish_tools::{MetaData, PafRecord, Summary, ReadfishSummary};
    ///
    /// let mut summary = Summary::new();
    /// let mut readfish_summary = ReadfishSummary::new(summary);
    ///
    /// // Create a MetaData object with relevant information
    /// let meta_data = MetaData {
    ///     condition_name: String::from("condition1"),
    ///     on_target: true,
    ///     paf_line: String::from("read123 200 0 200 + contig123 300 0 300 200 200 50 ch=1"),
    /// };
    ///
    /// // Update the summary using the MetaData
    /// let result = readfish_summary.update_summary(meta_data);
    /// assert!(result.is_ok());
    /// ```
    fn update_summary(&mut self, meta_data: MetaData) -> PyResult<()> {
        let paf_line = meta_data.paf_line;
        let t: Vec<&str> = paf_line.split_ascii_whitespace().collect();
        // Todo do without clone
        let paf_record = PafRecord::new(t).unwrap();
        {
            let mut x = self.summary.borrow_mut();
            let y = x.conditions(meta_data.condition_name.as_str());
            y.update(paf_record, meta_data.on_target).unwrap();
        }
        Ok(())
    }

    /// Add a target to teh condition and contig summary, summing up aggregated stats as we go
    pub fn add_target(
        &mut self,
        condition_name: String,
        contig: String,
        contig_len: usize,
        start: usize,
        end: usize,
    ) -> PyResult<()> {
        {
            let mut summary = self.summary.borrow_mut();
            let y = summary.conditions(condition_name.as_str());
            y.add_target(contig, contig_len, start, end).unwrap();
        }
        Ok(())
    }

    /// Set the total number of reads in the sequencing data.
    pub fn add_total_reads(
        &mut self,
        condition_name: String,
        contig_name: String,
        contig_len: usize,
        total_reads: usize,
    ) -> PyResult<()> {
        {
            let mut summary = self.summary.borrow_mut();
            let y = summary.conditions(condition_name.as_str());
            y.add_total_reads(total_reads);
            let contig_summary = y.get_or_add_contig(&contig_name, contig_len);
            contig_summary.add_total_reads(total_reads);
        }
        Ok(())
    }

    /// Prints the summary of the `ReadfishSummary` to the standard output.
    ///
    /// This method borrows the `ReadfishSummary` immutably and prints its summary to the standard output.
    /// The summary is obtained by calling the `borrow` method on the `RefCell<Summary>` field of the
    /// `ReadfishSummary`.
    ///
    /// # Returns
    ///
    /// This function returns a `PyResult<()>` to indicate success or failure. If the summary is
    /// successfully printed, `Ok(())` is returned. If an error occurs during printing, an appropriate
    /// `PyErr` will be set, and `Err` will be returned.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # use pyo3::prelude::*;
    /// # use std::cell::RefCell;
    /// # use std::collections::HashMap;
    /// # use std::path::PathBuf;
    /// # use std::error::Error;
    ///
    /// # #[pyclass]
    /// # pub struct ReadfishSummary {
    /// #     // Fields of ReadfishSummary
    /// #     // ...
    /// # }
    /// #
    /// # #[pymethods]
    /// # impl ReadfishSummary {
    /// #     #[getter]
    /// #     pub fn summary(&self) -> PyResult<Ref<Summary>> {
    /// #         unimplemented!()
    /// #     }
    /// #
    /// /// Method to print the summary of ReadfishSummary.
    /// pub fn print_summary(&self) -> PyResult<()> {
    ///     println!("{}", self.summary.borrow());
    ///     Ok(())
    /// }
    /// # }
    /// ```
    ///
    /// The method can be called on an instance of `ReadfishSummary` to print its summary.
    ///
    /// ```rust,ignore
    /// # use pyo3::prelude::*;
    /// # use std::cell::RefCell;
    /// # use std::collections::HashMap;
    /// # use std::path::PathBuf;
    /// # use std::error::Error;
    /// #
    /// # #[pyclass]
    /// # pub struct ReadfishSummary {
    /// #     // Fields of ReadfishSummary
    /// #     // ...
    /// # }
    /// #
    /// # #[pymethods]
    /// # impl ReadfishSummary {
    /// #     #[getter]
    /// #     pub fn summary(&self) -> PyResult<Ref<Summary>> {
    /// #         unimplemented!()
    /// #     }
    /// #
    /// #     /// Method to print the summary of ReadfishSummary.
    /// #     pub fn print_summary(&self) -> PyResult<()> {
    /// #         println!("{}", self.summary.borrow());
    /// #         Ok(())
    /// #     }
    /// # }
    /// #
    /// # fn main() -> PyResult<()> {
    /// #     Python::with_gil(|py| {
    /// #         let gil = Python::acquire_gil();
    /// #         let py = gil.python();
    /// #
    /// #         // Create an instance of ReadfishSummary and call the print_summary method
    /// #         let readfish_summary = ReadfishSummary { /* Initialize fields... */ };
    /// #         readfish_summary.print_summary()?;
    /// #         Ok(())
    /// #     })
    /// # }
    /// ```
    pub fn print_summary(&self) -> PyResult<()> {
        println!("{}", self.summary.borrow());
        Ok(())
    }
}
/// A Python module implemented in Rust.
#[pymodule]
fn readfish_summarise(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<ReadfishSummary>()?;
    m.add_class::<MetaData>()?;
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::path::PathBuf;
    fn get_resource_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources");
        path
    }

    fn get_test_file(file: &str) -> PathBuf {
        let mut path = get_resource_dir();
        path.push(file);
        path
    }

    #[test]
    fn test_update_lengths() {
        // Create a PAF record with a query length of 100
        let paf = PafRecord::new(
            "read123 100 0 100 + contig123 300 0 300 200 200 50 ch=1"
                .split(' ')
                .collect(),
        )
        .unwrap();

        // Create a MeanReadLengths instance
        let mut mean_lengths = MeanReadLengths::new();

        // Initially, all mean lengths should be zero
        assert_eq!(mean_lengths.on_target, 0);
        assert_eq!(mean_lengths.off_target, 0);
        assert_eq!(mean_lengths.total, 0);

        // Update with an on-target read
        mean_lengths.update_lengths(&paf, true);

        // After the update, only on_target and total should be updated
        assert_eq!(mean_lengths.on_target, 100);
        assert_eq!(mean_lengths.off_target, 0);
        assert_eq!(mean_lengths.total, 100);

        // Update with an off-target read
        mean_lengths.update_lengths(&paf, false);

        // After the update, off_target and total should be updated
        assert_eq!(mean_lengths.on_target, 100);
        assert_eq!(mean_lengths.off_target, 100);
        assert_eq!(mean_lengths.total, 100);
        // Create a PAF record with a query length of 100
        let paf = PafRecord::new(
            "read123 150 0 100 + contig123 300 0 300 200 200 50 ch=1"
                .split(' ')
                .collect(),
        )
        .unwrap();
        // Update with an off-target read with a different length
        mean_lengths.update_lengths(&paf, false);

        // After the update, off_target and total should be updated
        assert_eq!(mean_lengths.on_target, 100);
        assert_eq!(mean_lengths.off_target, 125);
        assert_eq!(mean_lengths.total, 116);
    }
}
