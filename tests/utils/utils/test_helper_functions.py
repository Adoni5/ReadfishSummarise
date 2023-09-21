"""
Write out a FASTQ file using a set randiom seed from the included truncated
yeast reference file that is only the first 8kn of each contig.
A stats.tsv file for each contig is written out and can be compared to the output
 of the command

`readfish stats --toml tests/static/stats_test/yeast_summary_test.toml
 --fastq-directory tests/static/stats_test/ --no-paf-out --no-demultiplex`
Where there is only one region in the toml file. The stats.tsv file and the
 simulated fastq file are written out to the current working directory.

"""
import random
from collections import defaultdict
from statistics import median

from pyfastx import Fasta

# Constants
REFERENCE_PATH = "yeast_ref_8kb_contigs.fna.gz"
NUM_READS = 10
MAX_READ_LENGTH = 8000
OUTPUT_FILE = "stats.tsv"
RANDOM_SEED = 42

# Setting the seed for reproducibility
random.seed(RANDOM_SEED)


def generate_reads(reference, num_reads, max_length):
    reads = defaultdict(list)
    for contig_name, contig_seq in Fasta(reference, build_index=False):
        contig_length = len(contig_seq)
        reads[contig_name].extend(
            [
                contig_seq[start : start + random.randint(500, max_length)]
                for start in random.sample(range(0, contig_length - 100), num_reads)
            ]
        )

    return reads


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    half_total = total / 2
    n50 = 0
    agg = 0
    for length in sorted_lengths:
        agg += length
        if agg >= half_total:
            n50 = length
            break
    return n50


def calculate_statistics(reads):
    stats = []
    print(reads)
    for contig, contig_reads in reads.items():
        lengths = [len(read) for read in contig_reads]
        yield_value = sum(lengths)
        median_length = median(lengths)
        n50_value = calculate_n50(lengths)
        stats.append((contig, median_length, yield_value, n50_value))
    return stats


def write_to_tsv(stats, output_file):
    with open(output_file, "w") as f:
        # Write header
        f.write("Contig\tMedian Length\tYield\tN50\n")
        for contig, median_length, yield_value, n50 in stats:
            f.write(f"{contig}\t{median_length}\t{yield_value}\t{n50}\n")


def write_fastq(reads, filename):
    with open(filename, "w") as f:
        for contig_name, contig_reads in reads.items():
            for i, read in enumerate(contig_reads, start=1):
                # Fake quality string of 'I' chars (very high quality)
                fake_quality = "I" * len(read)
                f.write(
                    f"@{contig_name}_simulated_read_{i} read={i}"
                    f" ch={random.randint(0, 512)}\n{read}\n+\n{fake_quality}\n"
                )


def main():
    reads = generate_reads(REFERENCE_PATH, NUM_READS, MAX_READ_LENGTH)
    stats = calculate_statistics(reads)
    write_to_tsv(stats, OUTPUT_FILE)
    write_fastq(reads, "simulated_reads.fastq")


if __name__ == "__main__":
    main()
