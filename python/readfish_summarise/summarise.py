"""
Stats for various files produced during readfish experiments.
"""
from itertools import chain
from pathlib import Path

import click
from readfish._config import Conf
from readfish.plugins.abc import AlignerABC
from readfish.plugins.mappy import UNMAPPED_PAF
from readfish.plugins.utils import Action, get_contig_lengths
from readfish_summarise.fastq_utils import (
    batched,
    is_fastq_file,
    update_summary,
    yield_reads_for_alignment,
)
from readfish_summarise.readfish_summarise import ReadfishSummary


@click.command()
@click.option(
    "--toml",
    help="The path to the experiment configuration TOML file.",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
)
@click.option(
    "--fastq-directory",
    help="The path to the experiment MinKnow data directory.",
    required=True,
    type=click.Path(exists=True, dir_okay=True, resolve_path=True, path_type=Path),
)
@click.option(
    "--demultiplex/--no-demultiplex",
    help="Switch demultiplex fastq reads by condition and unblocked/sequenced.",
    default=True,
)
@click.option(
    "--paf-out / --no-paf-out",
    help="Switch outputing a PAF file for each condition and unblocked/sequenced.",
    default=True,
)
@click.option(
    "--gzip/--no-gzip",
    help="Switch compressing demultiplexed FASTQ and PAF files. CURRENTLY UNUSED.",
    default=True,
)
def fastq(
    toml: str | Path,
    fastq_directory: str | Path | None = None,
    demultiplex: bool = True,
    paf_out: bool = True,
    gzip=True,
) -> None:
    """
    Create stats tables a previous Readfish experiment based on the
    readfish configuration TOML and FASTQ data from the experiment.

    :param toml: The path to the experiment configuration TOML file.
    :param fastq_directory: The path to the directory containing FASTQ files.
    :param demultiplex: Demultiplex the fastq files into
        `condition`_`sequenced/unblocked`.fastq, defaults to True
    :param paf_out: Write out alignments as they are created, defaults to True
    """
    _fastq(toml, fastq_directory, demultiplex, paf_out)


def _fastq(
    toml: str | Path,
    fastq_directory: str | Path | None = None,
    demultiplex: bool = True,
    paf_out: bool = True,
    *args,
    **kwargs,
):
    """
    Internal fastq function that can be called from other functions.
        Can also be imported.
        Create stats tables a previous Readfish experiment based on the
    readfish configuration TOML and FASTQ data from the experiment.

    :param toml:
    :param fastq_directory:


    :param toml: The path to the experiment configuration TOML file.
    :param fastq_directory: The path to the directory containing FASTQ files,
        defaults to None
    :param demultiplex: Demultiplex the fastq files into
        `condition`_`sequenced/unblocked`.fastq, defaults to True
    :param paf_out: Write out alignments as they are created, defaults to True
    """
    assert Path(fastq_directory).exists(), "Fastq directory does not exist"
    if fastq_directory.is_file():
        assert is_fastq_file(
            fastq_directory
        ), "Passed fastq directory is not a fastq file"
    # Stores fastq file writers if used
    fastq_files = {}
    conf = Conf.from_file(toml, channels=512)
    if paf_out:
        paf_writer = open("readfish_fastq_stats.paf", "w", buffering=8192)
    summary = ReadfishSummary()
    mapper: AlignerABC = conf.mapper_settings.load_object(
        "Aligner", readfish_config=conf
    )
    # todo - use utils function for genome length
    contig_lengths = get_contig_lengths(mapper.aligner)
    ref_len = sum(contig_lengths.values())
    # Used to work out if we are double targets that are the same on both
    # strands when we drop down to one strand
    double_stranded_targets = set()
    # Sort our conditions and get target values so we can do a
    # little fun summary of target coverage
    for condition in chain(conf.regions, conf.barcodes.values()):
        # Add all summary conditions
        summary.add_condition(condition.name, ref_len)
        # Add each contig for the range of targets
        for contig, contig_length in contig_lengths.items():
            summary.add_contig_to_condition(condition.name, contig, contig_length)

        for contig, targets in (
            (contig, targets)
            for contig_targets in condition.targets._targets.values()
            for contig, targets in contig_targets.items()
        ):
            # todo - switch to yield targets.
            for start, end in targets:
                if (start, end, contig) not in double_stranded_targets:
                    contig_len = contig_lengths[contig]
                    summary.add_target(
                        condition.name,
                        contig,
                        int(start),
                        int(end) if end != float("inf") else contig_len,
                        ref_len,
                    )
                    double_stranded_targets.add((start, end, contig))
        for action in Action:
            fastq_files[(condition.name, action.name)] = open(
                f"{condition.name}_{action.name}.fastq", "w", buffering=8192
            )

    for batch in batched(
        yield_reads_for_alignment(
            fastq_directory=fastq_directory,
        ),
        50000,
    ):
        for result in mapper.map_reads(iter(batch)):
            control, condition = conf.get_conditions(result.channel, result.barcode)

            # We don't get the channel regions (if they exist) if we also have
            #  barcodes, so we need a check in update summary to handle this
            region = conf.get_region(result.channel)
            # Action is correct however, even if we have regions and barcodes
            action = condition.get_action(result.decision)
            on_target = action.name == "stop_receiving" or control

            if demultiplex:
                fastq_files[(condition.name, action.name)].write(
                    str(result.basecall_data)
                )
            paf_line = f"{result.read_id}\t{len(result.seq)}\t{UNMAPPED_PAF}"
            # No map - so add that to the summary
            if not result.alignment_data:
                update_summary(result, summary, condition, region, on_target, paf_line)
            #             # Won't run without mapping data, so either the block
            # #above or this one will run
            for index, alignment in enumerate(result.alignment_data):
                if index == 0:
                    paf_line = f"{result.read_id}\t{len(result.seq)}\t{alignment}"
                    update_summary(
                        result,
                        summary,
                        condition,
                        region,
                        on_target,
                        paf_line,
                    )
            if paf_out:
                paf_writer.write(paf_line + "\n")
    # Write out the summary
    summary.print_summary()
    # Close all files
    if paf_out:
        paf_writer.close()
    for file in fastq_files.values():
        file.close()


if __name__ == "__main__":
    print("dead")