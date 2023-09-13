"""
Stats for various files produced during readfish experiments.
"""
from pathlib import Path

import click
from readfish._config import Conf
from readfish_summarise import ReadfishSummary


@click.command()
@click.option(
    "--toml",
    is_flag=True,
    help="The path to the experiment configuration TOML file.",
    required=True,
)
@click.option(
    "--toml",
    is_flag=True,
    help="The path to the experiment MinKnow data directory.",
    required=True,
)
def fastq(
    toml: str | Path,
    fastq_directory: str | Path | None = None,
) -> None:
    """
    Create stats tables a previous Readfish experiment based on the
    readfish configuration TOML and FASTQ data from the experiment.

    :param toml: The path to the experiment configuration TOML file.
    :param fastq_directory: The path to the directory containing FASTQ files.
                            If not provided (None), the default directory will be used.

    :returns: None

    :Example:

    .. code-block:: python

        self.aligner.summarise(fastq_directory="path/to/fastq_data")
    """
    assert Path(fastq_directory).exists(), "Fastq directory does not exist"
    Conf.from_file(toml)
    ReadfishSummary()
    # # todo - use utils function for genome length
    # contig_lengths = {seq: len(.aligner.seq(seq)) for seq in self.aligner.seq_names}
    # ref_len = sum(contig_lengths.values())
    # # Used to work out if we are double targets that are the same on both
    # # strands when we drop down to one strand
    # double_stranded_targets = set()
    # # Sort our conditions and get target values so we can do a
    # # little fun summary of target coverage
    # for condition in chain(self.config.regions, self.config.barcodes.values()):
    #     summary.add_condition(condition.name, ref_len)
    #     # Add each contig for the range of targets
    #     for contig, contig_length in contig_lengths.items():
    #         summary.add_contig_to_condition(condition.name, contig, contig_length)

    #     for contig, targets in (
    #         (contig, targets)
    #         for contig_targets in condition.targets._targets.values()
    #         for contig, targets in contig_targets.items()
    #     ):
    #         # todo - switch to yield targets.
    #         for start, end in targets:
    #             if (start, end, contig) not in double_stranded_targets:
    #                 contig_len = contig_lengths[contig]
    #                 summary.add_target(
    #                     condition.name,
    #                     contig,
    #                     int(start),
    #                     int(end) if end != float("inf") else contig_len,
    #                     ref_len,
    #                 )
    #                 double_stranded_targets.add((start, end, contig))


#     for batch in batched(
#         yield_reads_for_alignment(
#             fastq_directory=fastq_directory,
#         ),
#         50000,
#     ):
#         for result in self.map_reads(iter(batch)):
#             control, condition = self.config.get_conditions(
#                 result.channel, result.barcode
#             )

#             # We don't get the channel regions (if they exist) if we also have
#  barcodes, so we need this check
#             aregion = self.config.get_region(result.channel)
#             # Action is correct however, even if we have regions and barcodes
#             action = condition.get_action(result.decision)
#             on_target = action.name == "stop_receiving" or control
#             # No map - so add that to the summary
#             if not result.alignment_data:
#                 paf_line = f"{result.read_id}\t{len(result.seq)}\t{UNMAPPED_PAF}"
#                 m = MetaData(
#                     condition_name=condition.name,
#                     on_target=on_target,
#                     paf_line=paf_line,
#                 )
#                 summary.update_summary(m)
#                 if region and not isinstance(condition, Region):
#                     m.condition_name = region.name
#                     summary.update_summary(m)
#             # Won't run without mapping data, so either the block
# #above or this one will run
#             for index, alignment in enumerate(result.alignment_data):
#                 if index == 0:
#                     paf_line = f"{result.read_id}\t{len(result.seq)}\t{alignment}"
#                     m = MetaData(
#                         condition_name=condition.name,
#                         on_target=on_target,
#                         paf_line=paf_line,
#                     )
#                     summary.update_summary(m)
#                     # Check that are not duplicating the region, which would happen
# if we didn't have barcodes
#                     if region and not isinstance(condition, Region):
#                         m.condition_name = region.name
#                         summary.update_summary(m)

#     summary.print_summary()

# if __name__ = "__main__":
#     summarise()
