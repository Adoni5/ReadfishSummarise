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
use csv::WriterBuilder;
use itertools::Itertools;
use log::info;
use num_format::{Locale, ToFormattedString};
use prettytable::{color, row, Attr, Cell, Row, Table};
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use readfish_tools::nanopore::format_bases;
use readfish_tools::paf::PafRecord;
use readfish_tools::MeanReadLengths;
use std::{cell::RefCell, collections::HashMap, error::Error, ops::Deref};

/// Dynamic result type for holding either a generic value or an error
pub type DynResult<T> = Result<T, Box<dyn Error + 'static>>;

/// Colour of the on target columns in table
const FG_ON: Attr = Attr::ForegroundColor(color::BRIGHT_CYAN);
/// COlour of the off target columns in table
const FG_OFF: Attr = Attr::ForegroundColor(color::BRIGHT_YELLOW);
/// Colour of other columns
const FG_OTHER: Attr = Attr::ForegroundColor(color::BRIGHT_WHITE);

/// Calculates the N50 and median from a dataset of u32 numbers.
///
/// # Arguments
///
/// * `numbers` - A mutable reference to a vector of u32 numbers.
///
/// # Returns
///
/// A tuple containing the N50 `Option<u32>` and median `Option<f64>`.
/// If the dataset is empty, both values will be None.
///
/// # Examples
///
/// ```rust,ignore
/// let mut dataset = vec![100, 200, 300, 400, 500];
/// let (n50, median) = calculate_n50_median(&mut dataset);
///
/// assert_eq!(n50, Some(400));
/// assert_eq!(median, Some(300.0));
/// ```
fn calculate_n50_median(numbers: &mut Vec<u32>) -> (Option<u32>, Option<f64>) {
    if numbers.is_empty() {
        return (None, None);
    }

    numbers.sort(); // Sort in ascending order

    let total_sum: u32 = numbers.iter().sum();
    let half_sum = total_sum / 2;
    let n50 = numbers
        .iter()
        .scan(0, |current_sum, &number| {
            *current_sum += number;
            Some((*current_sum, number))
        })
        .find(|(current_sum, _number)| current_sum >= &half_sum)
        .unwrap()
        .1;

    let median = if numbers.len() % 2 == 1 {
        Some(numbers[numbers.len() / 2] as f64)
    } else {
        let mid = numbers.len() / 2;
        let median_val = (numbers[mid - 1] as f64 + numbers[mid] as f64) / 2.0;
        Some(median_val)
    };

    (Some(n50), median)
}
/// Divide too floats, returning none if the divisor is zero
fn safe_divide(dividend: f64, divisor: f64) -> Option<f64> {
    if divisor == 0.0 {
        None // Division by zero, return None or handle it as needed
    } else {
        Some(dividend / divisor)
    }
}

/// Can't include new lines in headers if we are writing a csv out
fn get_write_out_safe_headers(write_out: bool) -> (&'static str, &'static str, &'static str) {
    let mut print_strings = (
        "Median read\nlengths",
        "Number of\ntargets",
        "Estimated\ncoverage",
    );

    // Conditionally remove newlines based on the write_out variable
    if write_out {
        print_strings = (
            "Median read lengths",
            "Number of targets",
            "Estimated coverage",
        );
    }
    print_strings
}

// Define a common summary structure
/// A common struct for the contig and condition summaries, holding all shared fields
#[derive(Debug)]
struct BaseSummary {
    /// The name or identifier of the sequencing data.
    pub name: String,
    /// The total number of aligned reads in the sequencing data.
    pub mapped_reads: usize,
    /// The total number of unaligned reads in the sequencing data.
    pub unmapped_reads: usize,
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
    /// Number of target regions
    pub number_of_targets: usize,
    /// Number of bases covered by targets
    pub number_target_bases: usize,
    /// Length of the reference
    pub ref_length: usize,
    /// on target read lengths
    pub on_target_read_lengths: RefCell<Vec<u32>>,
    /// off target read lengths
    pub off_target_read_lengths: RefCell<Vec<u32>>,
}

impl BaseSummary {
    /// Return a new BaseSummary with default values for all fields except `name` and `ref_length`.
    fn new(name: String, ref_length: usize) -> Self {
        BaseSummary {
            name,
            mapped_reads: 0,
            unmapped_reads: 0,
            off_target_alignment_count: 0,
            on_target_alignment_count: 0,
            off_target_yield: 0,
            on_target_yield: 0,
            mean_read_lengths: MeanReadLengths::new(),
            number_of_targets: 0,
            number_target_bases: 0,
            ref_length,
            off_target_read_lengths: RefCell::new(Vec::new()),
            on_target_read_lengths: RefCell::new(Vec::new()),
        }
    }
}

/// A trait for summarizing contig or sequencing data.
///
/// This trait defines methods for summarizing various metrics related to contig or sequencing data,
/// including alignment counts, yield percentages, and mean read lengths. Implementing types should
/// either override or use the default implementations of these methods.
pub trait Summarise {
    /// Get the total number of alignments (on and off target combined).
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ContigSummary::new("Contig1".to_string(), 1000);
    /// assert_eq!(summary.total_alignments(), 42);
    /// ```
    fn total_alignments(&self) -> usize {
        self.off_target_alignment_count() + self.on_target_alignment_count()
    }

    /// Get the percentage of on-target alignments against total alignments.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.on_target_alignment_percent(), 84.0);
    /// ```
    fn on_target_alignment_percent(&self) -> f64 {
        self.on_target_alignment_count() as f64 / self.total_alignments() as f64 * 100.0
    }

    /// Get the percentage of off-target alignments against total alignments.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.off_target_alignment_percent(), 16.0);
    /// ```
    fn off_target_alignment_percent(&self) -> f64 {
        self.off_target_alignment_count() as f64 / self.total_alignments() as f64 * 100.0
    }
    /// Get the percentage of on-target yield against total yield.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.on_target_yield_percent(), 70.0);
    /// ```
    fn on_target_yield_percent(&self) -> f64 {
        self.on_target_yield() as f64 / self.total_yield() as f64 * 100.0
    }

    /// Get the percentage of off-target yield against total yield.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.off_target_yield_percent(), 30.0);
    /// ```
    fn off_target_yield_percent(&self) -> f64 {
        self.off_target_yield() as f64 / self.total_yield() as f64 * 100.0
    }
    /// Get the formatted on-target yield.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.on_target_yield_formatted(), "1.23 MB");
    /// ```
    fn on_target_yield_formatted(&self) -> String {
        format_bases(self.on_target_yield())
    }
    /// Get the formatted off-target yield.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.off_target_yield_formatted(), "456.78 KB");
    /// ```
    fn off_target_yield_formatted(&self) -> String {
        format_bases(self.off_target_yield())
    }
    /// Get the total yield (base pairs) of on and off-target reads.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.total_yield(), 123456);
    /// ```
    fn total_yield(&self) -> usize {
        self.on_target_yield() + self.off_target_yield()
    }
    /// Get the formatted total yield.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.total_yield_formatted(), "123.46 KB");
    /// ```
    fn total_yield_formatted(&self) -> String {
        format_bases(self.total_yield())
    }
    /// Get the yield ratio (on target:off target).
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.yield_ratio(), "1:3.14");
    /// ```
    fn yield_ratio(&self) -> String {
        format!(
            "{}:{:.2}",
            safe_divide(self.on_target_yield() as f64, self.on_target_yield() as f64)
                .unwrap_or(0.0),
            safe_divide(
                self.off_target_yield() as f64,
                self.on_target_yield() as f64
            )
            .unwrap_or(0.0)
        )
        .to_string()
    }
    /// Estimate on-target coverage.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.estimated_target_coverage(), "42.0% X");
    /// ```
    fn estimated_target_coverage(&self) -> String {
        let yieldy = safe_divide(
            self.on_target_yield() as f64,
            *self.number_target_bases() as f64,
        )
        .unwrap_or(0.0);
        format!("{:.2} X", yieldy)
    }

    /// Estimate the percentage of the total reference genome that is a target.
    ///
    /// # Example
    ///
    /// ```
    /// let summary = ConditionSummary::new("Condition1".to_string(), 5000);
    /// assert_eq!(summary.percent_of_genome_target(), "10.53%");
    /// ```
    fn percent_of_genome_target(&self) -> String {
        format!(
            "{:.2}%",
            safe_divide(*self.number_target_bases() as f64, self.length() as f64).unwrap_or(0.0)
                * 100.0
        )
        .to_string()
    }
    /// Mean read length of all reads on the contig.
    fn mean_read_length(&self) -> usize {
        self.mean_read_lengths().total as usize
    }
    /// On target mean read length of all reads on the contig.
    fn on_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths().on_target as usize
    }
    /// Off target mean read length of all reads on the contig.
    fn off_target_mean_read_length(&self) -> usize {
        self.mean_read_lengths().off_target as usize
    }

    // Define these as associated methods (require implementation in each struct)
    /// Add a readfish target to the Contig summary
    fn add_target(&mut self, start: usize, end: usize) -> DynResult<()>;

    /// Get the total number of reads
    fn total_reads(&self) -> usize;

    /// Add a count of mapped reads to the total count
    fn add_mapped_reads(&mut self, total_reads: usize);

    /// Add a count of unmapped reads to the total count
    fn add_unmapped_reads(&mut self, total_reads: usize);

    /// Get access to the mean read lengths tracker
    fn mean_read_lengths(&self) -> &MeanReadLengths;
    /// Get the off target alignment count
    fn off_target_alignment_count(&self) -> usize;
    /// Get the on target alignment count
    fn on_target_alignment_count(&self) -> usize;
    /// Get the on target yield
    fn on_target_yield(&self) -> usize;
    /// Get the off target yield
    fn off_target_yield(&self) -> usize;
    /// Get the number of target bases
    fn number_target_bases(&self) -> &usize;
    /// Get the length of the contig
    fn length(&self) -> usize;
    /// N50 and median on target read length
    fn on_target_n50_median(&self) -> (Option<u32>, Option<f64>);
    /// N50 and median off target read length
    fn off_target_n50_median(&self) -> (Option<u32>, Option<f64>);
    /// Get the combined n50 and median
    fn n50_median(&self) -> (Option<u32>, Option<f64>);
}

impl Summarise for BaseSummary {
    fn off_target_alignment_count(&self) -> usize {
        self.off_target_alignment_count
    }

    fn on_target_alignment_count(&self) -> usize {
        self.on_target_alignment_count
    }

    fn on_target_yield(&self) -> usize {
        self.on_target_yield
    }

    fn off_target_yield(&self) -> usize {
        self.off_target_yield
    }

    fn number_target_bases(&self) -> &usize {
        &self.number_target_bases
    }

    fn length(&self) -> usize {
        self.ref_length
    }

    fn mean_read_lengths(&self) -> &MeanReadLengths {
        &self.mean_read_lengths
    }

    fn total_reads(&self) -> usize {
        self.mapped_reads + self.unmapped_reads
    }

    fn add_target(&mut self, start: usize, end: usize) -> DynResult<()> {
        self.number_of_targets += 1;
        self.number_target_bases += end - start;
        Ok(())
    }

    fn add_mapped_reads(&mut self, total_reads: usize) {
        self.mapped_reads += total_reads;
    }

    fn add_unmapped_reads(&mut self, total_reads: usize) {
        self.unmapped_reads += total_reads;
    }

    fn on_target_n50_median(&self) -> (Option<u32>, Option<f64>) {
        let mut on_target_read_lengths = self.on_target_read_lengths.borrow_mut();
        calculate_n50_median(&mut on_target_read_lengths)
    }

    fn off_target_n50_median(&self) -> (Option<u32>, Option<f64>) {
        let mut off_target_read_lengths = self.off_target_read_lengths.borrow_mut();
        calculate_n50_median(&mut off_target_read_lengths)
    }

    fn n50_median(&self) -> (Option<u32>, Option<f64>) {
        let mut on_target_read_lengths = self.on_target_read_lengths.borrow_mut();
        let mut off_target_read_lengths = self.off_target_read_lengths.borrow_mut();
        off_target_read_lengths.append(&mut on_target_read_lengths);
        calculate_n50_median(&mut off_target_read_lengths)
    }
}
/// Represents a summary of a contig or sequence from a sequencing experiment.
/// It includes various metrics related to the contig's characteristics and read mapping.
#[derive(Debug)]
pub struct ContigSummary {
    /// Holds the base summary fields and methods
    data: BaseSummary,
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
            data: BaseSummary::new(name, length),
        }
    }
}
#[derive(Debug)]
/// Represents a summary of sequencing data, including various metrics related to the output of the experiment.
pub struct ConditionSummary {
    /// Holds the base summary fields and methods
    data: BaseSummary,
    /// A vector of `ContigSummary` representing summaries of individual contigs or sequences
    /// in the sequencing data.
    pub contigs: HashMap<String, ContigSummary>,
}

impl ConditionSummary {
    /// Add a target to the COnditions summary and the correct contig summary
    pub fn add_target(&mut self, contig: String, start: usize, end: usize) -> DynResult<()> {
        self.data.number_of_targets += 1;
        self.data.number_target_bases += end - start;
        let contig = self.get_contig(&contig).unwrap();
        contig.data.add_target(start, end)?;
        Ok(())
    }

    /// Add a contig to the COnditions summary
    pub fn add_contig(&mut self, contig: String, contig_len: usize) -> DynResult<()> {
        let _contig = self.get_or_add_contig(&contig, contig_len);
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
        self.data.mean_read_lengths.update_lengths(&paf, on_target);
        if on_target {
            self.data.on_target_alignment_count += 1;
            self.data.on_target_yield += paf.query_length;
            // self.on_target_mean_read_quality += paf.tlen as f64;
        } else {
            self.data.off_target_alignment_count += 1;
            self.data.off_target_yield += paf.query_length;
            // self.off_target_mean_read_quality += paf.tlen as f64;
        }

        // Check if it's equal to "*" - which is from the readfish unmapped PAF const.
        let target_name: String = if &paf.target_name == "*" {
            // If it's "*", create a new String
            self.data.add_unmapped_reads(1);
            String::from("unmapped")
        } else {
            // If it's not "*", use the original borrowed reference
            self.data.add_mapped_reads(1);
            paf.target_name.to_string() // Convert &str to String
        };

        let contig = self.get_or_add_contig(&target_name, paf.target_length)?;
        if contig.data.name == "unmapped" {
            contig.data.add_unmapped_reads(1);
        } else {
            contig.data.add_mapped_reads(1);
        }
        contig
            .data
            .mean_read_lengths
            .update_lengths(&paf, on_target);
        if on_target {
            contig.data.on_target_alignment_count += 1;
            contig.data.on_target_yield += paf.query_length;
            contig
                .data
                .on_target_read_lengths
                .borrow_mut()
                .push(paf.query_length as u32);
        } else {
            contig.data.off_target_alignment_count += 1;
            contig.data.off_target_yield += paf.query_length;
            contig
                .data
                .off_target_read_lengths
                .borrow_mut()
                .push(paf.query_length as u32);
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
    /// * `ref_length` - The length of the reference
    pub fn new(name: String, ref_length: usize) -> Self {
        ConditionSummary {
            data: BaseSummary::new(name, ref_length),
            contigs: HashMap::new(),
        }
    }

    /// Get a reference to the vector of `ContigSummary`.
    pub fn contigs(&self) -> &HashMap<String, ContigSummary> {
        &self.contigs
    }

    /// Get a mutable reference to the vector of `ContigSummary`.
    pub fn contigs_mut(&mut self) -> &mut HashMap<String, ContigSummary> {
        &mut self.contigs
    }

    /// Get mutable references to the contigs
    pub fn get_contig(&mut self, contig: &str) -> Option<&mut ContigSummary> {
        self.contigs.get_mut(contig)
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
    pub fn get_or_add_contig(
        &mut self,
        contig: &str,
        length: usize,
    ) -> DynResult<&mut ContigSummary> {
        Ok(self
            .contigs
            .entry(contig.to_string())
            .or_insert(ContigSummary::new(contig.to_string(), length)))
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

impl Summary {
    /// Create a new `Summary` instance with default values for all fields.
    fn new() -> Self {
        Summary {
            conditions: HashMap::new(),
        }
    }
    ///ahhh
    fn display(&self) -> DynResult<()> {
        let condition_table = self.create_condition_table(false);
        condition_table.printstd();
        println!("Contigs:");

        for condition_summary in self
            .conditions
            .values()
            .sorted_by(|key1, key2| natord::compare(&key1.data.name, &key2.data.name))
        {
            let mut contig_table = Table::new();
            contig_table = self.create_contig_table(contig_table, condition_summary, false, 0);
            contig_table.printstd();
        }
        Ok(())
    }
    ///aaa
    fn log_table(&self) -> DynResult<()> {
        let condition_table = self.create_condition_table(false);
        let mut tables: Vec<String> = vec![];
        tables.push(condition_table.to_string());
        for condition_summary in self
            .conditions
            .values()
            .sorted_by(|key1, key2| natord::compare(&key1.data.name, &key2.data.name))
        {
            let mut contig_table = Table::new();
            contig_table = self.create_contig_table(contig_table, condition_summary, false, 0);
            tables.push(contig_table.to_string());
        }
        info!("{}", tables.join("\n\n"));
        Ok(())
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
        ref_length: usize,
    ) -> &mut ConditionSummary {
        self.conditions
            .entry(condition_name.to_string())
            .or_insert(ConditionSummary::new(
                condition_name.to_string(),
                ref_length,
            ))
    }

    /// Get a mutable reference to the `ConditionSummary` associated with the specified condition name.
    ///
    /// If a `ConditionSummary` with the provided `condition_name` exists in the summary,
    /// this method returns a mutable reference to it. If no matching `ConditionSummary`
    /// is found, it will panic with an error message.
    ///
    /// # Arguments
    ///
    /// * `condition_name` - The name or identifier of the condition to retrieve.
    ///
    /// # Panics
    ///
    /// This method will panic if the specified `condition_name` is not found in the summary.
    ///
    /// # Example
    ///
    /// ```rust
    /// use your_crate::Summary;
    ///
    /// let mut summary = Summary::new();
    /// let condition_name = "ConditionA";
    /// let condition_summary = summary.get_condition(condition_name);
    ///
    /// // Now you can modify the `ConditionSummary` fields or access its properties
    /// println!("Condition name: {}", condition_summary.name);
    /// ```
    pub fn get_condition<T: Deref<Target = str> + std::fmt::Debug>(
        &mut self,
        condition_name: T,
    ) -> &mut ConditionSummary {
        self.conditions
            .get_mut(&condition_name.to_string())
            .unwrap_or_else(|| panic!("Unable to find Condition in summary: {:#?}", condition_name))
    }
    /// Write out table to svc file
    pub fn to_csv(&self, file_path: impl AsRef<str>) -> DynResult<()> {
        let condition_table = self.create_condition_table(true);
        let mut contig_table = Table::new();
        for (index, condition_summary) in self
            .conditions
            .values()
            .sorted_by(|key1, key2| natord::compare(&key1.data.name, &key2.data.name))
            .enumerate()
        {
            contig_table = self.create_contig_table(contig_table, condition_summary, true, index);
        }
        let wtr = WriterBuilder::new()
            .flexible(true)
            .from_path(format!("{}_conditions.csv", file_path.as_ref()))?;
        condition_table.to_csv_writer(wtr)?;
        let wtr = WriterBuilder::new()
            .flexible(true)
            .from_path(format!("{}_contigs.csv", file_path.as_ref()))?;
        contig_table.to_csv_writer(wtr)?;
        Ok(())
    }

    /// Create the contig table - takes in the table and the condition summary
    pub fn create_contig_table(
        &self,
        mut contig_table: Table,
        condition_summary: &ConditionSummary,
        write_out: bool,
        index: usize,
    ) -> Table {
        let print_strings = get_write_out_safe_headers(write_out);
        if !write_out {
            contig_table.add_row(Row::new(vec![
                Cell::new("Condition Name")
                    .with_style(Attr::Bold)
                    .with_style(FG_OTHER),
                Cell::new(&condition_summary.data.name)
                    .with_style(Attr::Standout(true))
                    .with_style(FG_OTHER)
                    .with_hspan(21),
            ]));
        }
        let mut header_cells = vec![
            Cell::new("Condition")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Contig")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Contig Length")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Reads")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Alignments")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("Yield")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(4),
            Cell::new(print_strings.0)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new("N50")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new(print_strings.1)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Percent target")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new(print_strings.2)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
        ];

        // duplicate the cells so if we are writing out we write out valid CSV
        // Duplicate elements as many times as their hspan
        let empty_cell_filler = if write_out { "Total" } else { "" };
        if write_out {
            let index_to_duplicate = [3, 4, 5, 6, 7];
            for index in index_to_duplicate.iter().rev() {
                let cell_to_duplicate = header_cells[*index].clone();
                for _ in 1..cell_to_duplicate.get_hspan() {
                    header_cells.insert(*index, cell_to_duplicate.clone().with_hspan(1));
                }
            }
        }
        // If we are writing out one big table we don't want to rewrite the header line each condition
        // so only write it out if it is the first iteration.
        if index == 0 {
            contig_table.add_row(Row::new(header_cells));
            contig_table.add_row(row!["", "", empty_cell_filler, FBb->"Mapped", FYb->"Unmapped", b->"Total", FBb->"On-Target", FYb->"Off-Target", b->"Total",  FBb->"On-Target", FYb->"Off-Target", b->"Total", b->"Ratio", FBb->"On-target", FYb->"Off-target", b->"Combined", FBb->"On-Target", FYb->"Off-Target", b->"Total",empty_cell_filler, empty_cell_filler, empty_cell_filler]);
        }
        for (contig_name, contig_summary) in condition_summary
            .contigs
            .iter()
            .sorted_by(|(key1, _), (key2, _)| natord::compare(key1, key2))
        {
            let (on_target_n50, on_target_median) = contig_summary.data.on_target_n50_median();
            let (off_target_n50, off_target_median) = contig_summary.data.off_target_n50_median();
            let (n50, median) = contig_summary.data.n50_median();
            let cells = vec![
                // Condition name
                Cell::new(&condition_summary.data.name)
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
                // contig name
                Cell::new(contig_name)
                    .with_style(Attr::Blink)
                    .with_style(FG_OTHER),
                // contig length
                Cell::new(
                    &contig_summary
                        .data
                        .length()
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OTHER),
                // Number of mapped reads
                Cell::new(
                    &contig_summary
                        .data
                        .mapped_reads
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_ON),
                // Number of unmapped reads
                Cell::new(
                    &contig_summary
                        .data
                        .unmapped_reads
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OFF),
                // Number of reads
                Cell::new(
                    &contig_summary
                        .data
                        .total_reads()
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OTHER),
                // on target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    contig_summary
                        .data
                        .on_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    contig_summary.data.on_target_alignment_percent()
                ))
                .with_style(FG_ON),
                // off target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    contig_summary
                        .data
                        .off_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    contig_summary.data.off_target_alignment_percent()
                ))
                .with_style(FG_OFF),
                // total alignments
                Cell::new(
                    &contig_summary
                        .data
                        .total_alignments()
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OTHER),
                // on target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    contig_summary.data.on_target_yield_formatted(),
                    contig_summary.data.on_target_yield_percent()
                ))
                .with_style(FG_ON),
                //off target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    contig_summary.data.off_target_yield_formatted(),
                    contig_summary.data.off_target_yield_percent()
                ))
                .with_style(FG_OFF),
                // total yield
                Cell::new(&contig_summary.data.total_yield_formatted()).with_style(FG_OTHER),
                // yield ratio
                Cell::new(&contig_summary.data.yield_ratio()).with_style(FG_OTHER),
                // on target median read length
                Cell::new(&format_bases(on_target_median.unwrap_or(0_f64) as usize))
                    .with_style(FG_ON),
                // off target median read length
                Cell::new(&format_bases(off_target_median.unwrap_or(0_f64) as usize))
                    .with_style(FG_OFF),
                // median read length
                Cell::new(&format_bases(median.unwrap_or(0_f64) as usize)).with_style(FG_OTHER),
                // on target median read length
                Cell::new(&format_bases(on_target_n50.unwrap_or(0) as usize)).with_style(FG_ON),
                // off target median read length
                Cell::new(&format_bases(off_target_n50.unwrap_or(0) as usize)).with_style(FG_OFF),
                // median read length
                Cell::new(&format_bases(n50.unwrap_or(0) as usize)).with_style(FG_OTHER),
                // number of targets
                Cell::new(
                    &contig_summary
                        .data
                        .number_of_targets
                        .to_formatted_string(&Locale::en)
                        .to_string(),
                )
                .with_style(FG_OTHER),
                // Percent of contig that is a target
                Cell::new(&contig_summary.data.percent_of_genome_target()).with_style(FG_OTHER),
                // estimated target coverage
                Cell::new(&contig_summary.data.estimated_target_coverage()).with_style(FG_OTHER),
            ];

            contig_table.add_row(Row::new(cells));
        }
        contig_table
    }

    /// Create the condition table for printing / writing out
    fn create_condition_table(&self, write_out: bool) -> Table {
        // Todo rewrite to use Macro!
        let mut condition_table = Table::new();
        let print_strings = get_write_out_safe_headers(write_out);
        let mut header_cells = vec![
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
            Cell::new(print_strings.0)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN))
                .with_hspan(3),
            Cell::new(print_strings.1)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new("Percent target")
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
            Cell::new(print_strings.2)
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::BRIGHT_GREEN)),
        ];

        // duplicate the cells so if we are writing out we write out valid CSV
        // Duplicate elements as many times as their hspan
        let empty_cell_filler = if write_out { "Total" } else { "" };
        if write_out {
            let index_to_duplicate = [2, 3, 4];
            for index in index_to_duplicate.iter().rev() {
                let cell_to_duplicate = header_cells[*index].clone();
                for _ in 1..cell_to_duplicate.get_hspan() {
                    header_cells.insert(*index, cell_to_duplicate.clone().with_hspan(1));
                }
            }
        }
        let sub_header_row = row!["", empty_cell_filler, FBb->"On-Target", FYb->"Off-Target", FWb->"Total",  FBb->"On-Target", FYb->"Off-Target", FWb->"Total", FWb->"Ratio", FBb->"On-target", FYb->"Off-target", FWb->"Combined", empty_cell_filler, empty_cell_filler];
        condition_table.add_row(Row::new(header_cells.clone()));
        condition_table.add_row(sub_header_row.clone());
        for (condition_name, condition_summary) in self
            .conditions
            .iter()
            .sorted_by(|(key1, _), (key2, _)| natord::compare(key1, key2))
        {
            condition_table.add_row(Row::new(vec![
                Cell::new(condition_name).with_style(Attr::ForegroundColor(color::BRIGHT_YELLOW)),
                // total reads
                Cell::new(
                    &condition_summary
                        .data
                        .total_reads()
                        .to_formatted_string(&Locale::en),
                )
                .with_style(FG_OTHER),
                // on target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .data
                        .on_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    condition_summary.data.on_target_alignment_percent()
                ))
                .with_style(FG_ON),
                // off target alignment
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .data
                        .off_target_alignment_count
                        .to_formatted_string(&Locale::en),
                    condition_summary.data.off_target_alignment_percent()
                ))
                .with_style(FG_OFF),
                // total alignments
                Cell::new(
                    &condition_summary
                        .data
                        .total_alignments()
                        .to_formatted_string(&Locale::en)
                        .to_string(),
                )
                .with_style(FG_OTHER),
                // on target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary.data.on_target_yield_formatted(),
                    condition_summary.data.on_target_yield_percent()
                ))
                .with_style(FG_ON),
                // off target yield
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary.data.off_target_yield_formatted(),
                    condition_summary.data.off_target_yield_percent()
                ))
                .with_style(FG_OFF),
                // total yield
                Cell::new(&condition_summary.data.total_yield_formatted()).with_style(FG_OTHER),
                Cell::new(&condition_summary.data.yield_ratio()).with_style(FG_OTHER),
                // on target mean read length
                Cell::new(&format_bases(
                    condition_summary.data.on_target_mean_read_length(),
                ))
                .with_style(FG_ON),
                // off target mean read length
                Cell::new(&format_bases(
                    condition_summary.data.off_target_mean_read_length(),
                ))
                .with_style(FG_OFF),
                // mean read length
                Cell::new(&format_bases(condition_summary.data.mean_read_length()))
                    .with_style(FG_OTHER),
                // number of targets
                Cell::new(
                    &condition_summary
                        .data
                        .number_of_targets
                        .to_formatted_string(&Locale::en)
                        .to_string(),
                )
                .with_style(FG_OTHER),
                // Percent of genome that is a target
                Cell::new(&condition_summary.data.percent_of_genome_target()).with_style(FG_OTHER),
                // Estimated target coverage
                Cell::new(&condition_summary.data.estimated_target_coverage()).with_style(FG_OTHER),
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
        if !write_out {
            condition_table.add_row(sub_header_row);
            condition_table.add_row(Row::new(header_cells));
        };
        condition_table
    }
}

// PYTHON PyO3 STuff below ////////////////////////

#[pyclass]
#[derive(Debug, Clone)]
/// Stores metadata about a read's mapping and condition.
pub struct MetaData {
    /// The name of the condition from ReadFish analysis.
    #[pyo3(get, set)]
    pub condition_name: String,
    /// Indicates whether the read mapped to an on-target region of the genome.
    #[pyo3(get, set)]
    pub on_target: bool,
    /// The Pafline to be analysed
    #[pyo3(get, set)]
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

    /// Repr for the class
    fn __repr__(&self) -> String {
        // We use the `format!` macro to create a string. Its first argument is a
        // format string, followed by any number of parameters which replace the
        // `{}`'s in the format string.
        //
        format!(
            "Metadata({}, {}, {})",
            self.condition_name, self.on_target, self.paf_line
        )
    }

    /// String for the class
    fn __str__(&self) -> String {
        format!(
            "{0}, {1}, {2}",
            self.condition_name, self.on_target, self.paf_line
        )
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
    /// * `ref_length` - The length of the reference genome.
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
    #[pyo3(signature = (meta_data))]
    fn update_summary(&mut self, meta_data: MetaData) -> PyResult<()> {
        let paf_line = meta_data.paf_line;
        let t: Vec<&str> = paf_line.split_ascii_whitespace().collect();
        // Todo do without clone
        let paf_record = PafRecord::new(t).unwrap();
        {
            let mut x = self.summary.borrow_mut();
            let y = x.get_condition(meta_data.condition_name.as_str());
            y.update(paf_record, meta_data.on_target).unwrap();
        }
        Ok(())
    }
    /// Add a new condition to the summary.
    ///
    /// This method creates a new `ConditionSummary` with the provided `condition_name` and `ref_length`,
    /// and adds it to the summary. If a condition with the same name already exists, it will be replaced.
    ///
    /// # Arguments
    ///
    /// * `condition_name` - The name or identifier of the condition to add.
    /// * `ref_length` - The length of the reference associated with the condition.
    ///
    /// # Returns
    ///
    /// This function returns a `PyResult` indicating the success or failure of the operation.
    /// If the operation is successful, it returns an `Ok(())`. Otherwise, it returns an `Err` with an error message.
    ///
    /// # Example
    ///
    /// ```rust
    /// use your_crate::Summary;
    ///
    /// let mut summary = Summary::new();
    /// let condition_name = "ConditionA".to_string();
    /// let ref_length = 1000;
    ///
    /// // Add a new condition to the summary
    /// summary.add_condition(condition_name, ref_length).expect("Failed to add condition");
    /// ```
    pub fn add_condition(&mut self, condition_name: String, ref_length: usize) -> PyResult<()> {
        {
            let mut summary = self.summary.borrow_mut();
            summary.conditions(condition_name.as_str(), ref_length);
        }
        Ok(())
    }

    /// Add a contig to the contig summary, summing up aggregated stats as we go
    pub fn add_contig_to_condition(
        &mut self,
        condition_name: String,
        contig: String,
        contig_len: usize,
    ) -> PyResult<()> {
        {
            let mut summary = self.summary.borrow_mut();
            let y = summary.get_condition(condition_name.as_str());
            y.add_contig(contig, contig_len).unwrap();
        }
        Ok(())
    }

    /// Add a target to teh condition and contig summary, summing up aggregated stats as we go
    pub fn add_target(
        &mut self,
        condition_name: String,
        contig: String,
        start: usize,
        end: usize,
        ref_length: usize,
    ) -> PyResult<()> {
        {
            let mut summary = self.summary.borrow_mut();
            let y = summary.conditions(condition_name.as_str(), ref_length);
            y.add_target(contig, start, end).unwrap();
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
    #[pyo3(signature = (write_out=true))]
    pub fn print_summary(&self, write_out: bool) -> PyResult<()> {
        self.summary.borrow().display().unwrap();
        if write_out {
            let summary = self.summary.borrow();
            summary.to_csv("test.csv").unwrap();
        }
        Ok(())
    }

    /// Retrieves the summary as a string for printing  in python, loses colours
    pub fn log_summary(&self) -> PyResult<()> {
        self.summary.borrow().log_table().unwrap();
        Ok(())
    }
}
/// A Python module implemented in Rust.
#[pymodule]
fn readfish_summarise(_py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();
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

    #[cfg(test)]
    #[test]
    fn test_calculate_n50_median() {
        // Test case 1: Empty dataset
        let mut empty_dataset: Vec<u32> = vec![];
        let (n50_empty, median_empty) = calculate_n50_median(&mut empty_dataset);
        assert_eq!(n50_empty, None);
        assert_eq!(median_empty, None);

        // Test case 2: Odd-length dataset
        let mut odd_length_dataset: Vec<u32> = vec![10, 20, 30, 40, 50];
        let (n50_odd, median_odd) = calculate_n50_median(&mut odd_length_dataset);
        assert_eq!(n50_odd, Some(40));
        assert_eq!(median_odd, Some(30.0));

        // Test case 3: Even-length dataset
        let mut even_length_dataset: Vec<u32> = vec![10, 20, 30, 40];
        let (n50_even, median_even) = calculate_n50_median(&mut even_length_dataset);
        assert_eq!(n50_even, Some(30));
        assert_eq!(median_even, Some(25.0));

        // Test case 4: Dataset with duplicate values
        let mut dataset_with_duplicates: Vec<u32> = vec![10, 10, 20, 30, 30, 30, 40, 40, 50];
        let (n50_duplicates, median_duplicates) =
            calculate_n50_median(&mut dataset_with_duplicates);
        assert_eq!(n50_duplicates, Some(30));
        assert_eq!(median_duplicates, Some(30.0));

        // Test case 5: Not evenly spaced dataset
        let mut odd_length_dataset: Vec<u32> = vec![1, 78, 12, 3, 108, 1076];
        let (n50_random, median_uneven) = calculate_n50_median(&mut odd_length_dataset);
        assert_eq!(n50_random, Some(1076));
        assert_eq!(median_uneven, Some(45.0));

        // Test case 6: Dataset with unordered values
        let mut random_dataset: Vec<u32> = vec![30, 10, 29, 3, 7, 10000, 4, 23, 1];
        let (n50_random, median_random) = calculate_n50_median(&mut random_dataset);
        assert_eq!(n50_random, Some(10000));
        assert_eq!(median_random, Some(10.0));
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
