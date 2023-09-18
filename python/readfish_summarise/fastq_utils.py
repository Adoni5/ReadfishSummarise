from __future__ import annotations

import re
from io import TextIOWrapper
from itertools import islice
from pathlib import Path
from typing import Iterable, NamedTuple

from mappy import fastx_read
from readfish._config import Action, Barcode, Region
from readfish.plugins.utils import Result
from readfish_summarise.readfish_summarise import MetaData, ReadfishSummary


class FastqRecord(NamedTuple):
    name: str
    description: str
    sequence: str
    quality: str
    comment: str = "+"

    def __str__(self):
        fastq_string = "\n".join(
            [
                f"@{self.name} {self.description}",
                self.sequence,
                self.comment,
                self.quality,
            ]
        )
        return f"{fastq_string}\n"


def batched(iterable, n):
    """Batch data into tuples of length n. The last batch may be shorter."""
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def is_fastq_file(file_path: Path) -> bool:
    """
    Check a file suffix indicates a fastq file
    """
    types = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
    return bool(
        set(["".join(list(map(str.lower, file_path.suffixes)))]).intersection(types)
    )


def get_fq(directory: str | Path):
    """
    Given a directory, return a generator of fastq files.

    Parameters:

    :param directory (str or Path): The directory path to search for fastq files.

    Yields:

    :yield str: A path to a fastq file found in the given directory
        or its subdirectories.

    Example:
    --------

    .. code-block:: python

        for file_path in get_fq("resouces"):
            print(file_path)

    Output
    ------
    .. code-block:: python

        /path/to/directory/sample1.fastq
        /path/to/directory/sample2.fastq.gz

    Note:
        The function searches for files with extensions .fastq, .fastq.gz, .fq, .fq.gz
        in the specified directory and its subdirectories.
    """
    files = (str(p.resolve()) for p in Path(directory).glob("**/*") if is_fastq_file(p))
    yield from files


def yield_reads_for_alignment(fastq_directory: str | Path) -> Iterable[Result]:
    """
    Yield reads for alignment.

    This function yields reads for alignment from a specified fastq directory.

    :param fastq_directory: The path to the fastq directory.
    :return: An iterable of Result objects for alignment.

    :Example:

    .. code-block:: python

       # Assuming valid inputs and imports
       al = Alignment()  # Your alignment object

       # Iterate through reads from the fastq directory
       for read_data in yield_reads_for_alignment("/path/to/fastq/directory"):
           channel = read_data.channel
           read_number = read_data.read_number
           read_id = read_data.read_id
           seq = read_data.seq
           barcode = read_data.barcode
           # Perform alignment or other processing on the read data
           ...
    """
    # Define a regex pattern to capture each side of the "=" sign
    pattern = r"(\w+)=([^ =]+)"
    pattern = re.compile(pattern)
    for file in get_fq(fastq_directory):
        for name, seq, qual, comment in fastx_read(file, read_comment=True):
            # Find all matches of the pattern in the input string
            comments = dict(pattern.findall(comment))
            channel = int(comments["ch"])
            read_number = int(comments["read"])
            barcode = comments.get("barcode", None)
            yield Result(
                channel=channel,
                read_number=read_number,
                read_id=name,
                seq=seq,
                barcode=barcode,
                basecall_data=FastqRecord(
                    name=name, description=comment, sequence=seq, quality=qual
                ),  # res,
            )


def update_summary(
    result: Result,
    summary: ReadfishSummary,
    condition: Barcode | Region,
    region: Barcode | Region,
    on_target: bool,
    paf_line: str,
) -> bool:
    """
    Update the summary information for a given FASTQ read result.

    This function updates the provided summary with metadata regarding
    the alignment and target condition of a read.

    :param result: The FASTQ read result containing read details.
    :param summary: The summary object to be updated.
    :param condition: The condition for which the read was checked.
    :param region: The specific genomic region of interest for the read.
    :param on_target: Flag indicating if the read was on target.
    :param paf_line: The alignment paf line for the read.

    :note: If a region is provided and the condition is not of type 'Region',
           the function updates the summary with the barcode and also
           again with the region.

    :return: True if the summary was updated with a Barcode AND a region,
             False otherwise.
    """
    m = MetaData(
        condition_name=condition.name,
        on_target=on_target,
        paf_line=paf_line,
    )
    summary.update_summary(m)
    # Check that are not duplicating the region, which would happen
    # if we didn't have barcodes
    if region and not isinstance(condition, Region):
        m.condition_name = region.name
        summary.update_summary(m)
        return True
    return False


def write_out_fastq(
    control: bool,
    condition: Barcode | Region,
    action: Action,
    result: Result,
    fastq_files: dict[(str, str), TextIOWrapper],
):
    """
    Writes out FASTQ data based on the condition and action parameters.

    This function takes in the condition and action for a particular read and
    writes out the FASTQ data to the appropriate file. If the read is a control
    read without any targets, the function ensures that the action is labeled
    'stop_receiving' regardless of the input action.

    :param control: A flag indicating whether the read is from a channel in a
      control region
    :type control: bool

    :param condition: Specifies the condition (either Barcode or Region)
      under which the read falls.
    :type condition: Barcode | Region

    :param action: The determined action for the read, such as 'stop_receiving',
      'unblock', etc.
    :type action: Action

    :param result: Contains details about the read including its basecall data.
    :type result: Result

    :param fastq_files: A dictionary mapping conditions and actions to their respective
                        file output streams.
    :type fastq_files: dict[(str, str), TextIOWrapper]

    :return: None
    """
    # Control with no targets always gives an unblock decision,
    # which is incorrect so label stop receiving
    if control:
        action = Action.stop_receiving
    fastq_files[(condition.name, action.name)].write(str(result.basecall_data))
