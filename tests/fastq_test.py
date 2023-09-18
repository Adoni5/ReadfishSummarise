import re
from pathlib import Path

import pytest
from click.testing import CliRunner
from readfish.plugins.utils import Decision, Result
from readfish_summarise.fastq_utils import (
    FastqRecord,
    get_fq,
    yield_reads_for_alignment,
)
from readfish_summarise.summarise import fastq


def remove_ansi_escape_sequences(text):
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


TEST_FILES = Path(__file__).parent.resolve() / "static" / "summarise_test"

MMI_FILE = TEST_FILES / "yeast_8kb_contigs_test.mmi"
NON_BARCODED_FASTQ_DIRECTORY = TEST_FILES / "non_barcoded"
BARCODED_FASTQ_DIRECTORY = TEST_FILES / "barcoded"
LARGE_FILE_DIRECTORY = TEST_FILES / "large_file"


@pytest.fixture
def fastq_directory():
    return TEST_FILES / "non_barcoded"


@pytest.fixture
def barcoded_fastq_directory():
    return TEST_FILES / "barcoded"


@pytest.fixture
def mmi_file():
    return str(MMI_FILE)


def _generate_test_params():
    yield from [
        (
            TEST_FILES / "yeast_summary_test.toml",
            NON_BARCODED_FASTQ_DIRECTORY,
            TEST_FILES / "expected_summary.txt",
        ),
        (
            TEST_FILES / "yeast_summary_barcoded_test.toml",
            BARCODED_FASTQ_DIRECTORY,
            TEST_FILES / "expected_barcoded_summary.txt",
        ),
        (
            TEST_FILES / "yeast_summary_barcoded_regions_test.toml",
            BARCODED_FASTQ_DIRECTORY,
            TEST_FILES / "expected_barcoded_regions_summary.txt",
        ),
        (
            TEST_FILES / "yeast_summary_large_file_regions_test.toml",
            LARGE_FILE_DIRECTORY,
            TEST_FILES / "expected_summary_large_file.txt",
        ),
    ]


@pytest.mark.parametrize("toml,read_directory,expected", _generate_test_params())
def test_summarise(capfd, toml, read_directory, expected):
    runner = CliRunner()
    runner.invoke(
        fastq,
        [
            "--toml",
            str(toml),
            "--fastq-directory",
            str(read_directory),
            "--no-demultiplex",
            "--no-paf-out",
        ],
    )
    out, _err = capfd.readouterr()
    x = remove_ansi_escape_sequences(out)
    with open(expected, "rt") as fh:
        expected_message = fh.read()
    assert expected_message == x


def test_yield_reads_for_alignment(fastq_directory):
    alignments = 0

    alignments = sorted(
        list(yield_reads_for_alignment(fastq_directory)), key=lambda x: x.read_id
    )
    num_alignments = len(alignments)
    first_alignment = alignments[0]
    s = (
        "TAATATATATTTAATATAATTAAATATTATATTATATTATATTATATTATTTATTAAAAAAAAATCTATTACT"
        "TATTTTTTTTATTAATATATAAATTATTTATATAATTTATCATTTTTATTTATATATTATTATTTTTTATATATAAA"
        "TTAATATATATATATATTATATATACTTTTTTTTTTATAATATATCTATATATATAAATAAATATATTATATTATAT"
        "TTTTATATAATATATTATTAATTATTATTTTAATTTTCTATTCTATTGTGGGGGTCCCAATTATTATTTTCAATAAT"
        "AATTATTATTGGGACCCGGATATCTTCTTGTTTATCATTTATTATTTTATTAAATTTATTATTATTTTTAATTTATA"
        "TTTATATTATATAATTAATTATATCGTTTATACTCCTTCGGGGTCCCCGCCGGGGCGGGGACTTTATATTTTATTAT"
        "ATAATATATTATATTCTTATAATATATTTATTGATTATGTTATAAAATTTATTCTATGTGTGCTCTATATATATTTA"
        "ATATTCTGGTTATTATCACCCACCCCCTCCCCCTATTACGTCTCCGAGGTCCCGGTTTCGTAAGAAACCGGGACTTA"
        "TATATTTATAAATATAAATCTAACTTAATTAATAATTTAAATAATATACTTTATATTTTATAAATAAAAATAATTAT"
        "AACCTTTTTTATAATTATATATAATAATAATATATATTATCAAATAATTATTATTTCTTTTTTTTCTTTAATTAATT"
        "AATTAATTAATATTTTATAAAAATATATTTCTCCTTACGGGGTTCCGGCTCCCGTAGCCGGGGCCCGAAACTAAATA"
        "AAATATATTATTAATAATATTATATAATATAATAATAATATAATAATTTTATATAAATATATATTTATATATTAAAT"
        "TAAATTATAATTTTATTATGAAAATTATATCTTTTTTTTATATTTTTATATAATAAAAATATGTTATATATATATTA"
        "ATAATAAAAGGTAGTGAGGATTAAATAAATTATATAATAATTATAACT"
    )
    q = (
        "V2?t2`kB*''K;yTF%5$}M`ZY$~+*qmkDN|~Z+=y.Oiwvzr\\|M-.WW;A!Nh}!wm:GGdB"
        "]{K$!<YxqTlYS1~V,=7FYr-Nu[^}x<Op?5?Npo\\f\\pIt@~tDgm/uL~9&E'2.Ak>arc"
        "E1@YH:T[n;_?XcWwq`v{]FX@N%X%MU]'j9|Mg:e:6!~P-AZclXzR.V?j>8?jj6/CyE>A."
        'nu=mtpl^Mf\\)2ctyni"psmT!17#\\i0s&S|{K]sB>}jw\\PGJs%nugWG)YJc(@C{>)j?U'
        "NI`%dM8;!+'}qdZ+n!ltOOH]TR<Tt22d48WH#y@p-zSvm\\MtuF^&_7*UOxfFx`j5'@&]"
        "\\G3V7QSBy%}pQqb1SV4t!GyC?*.{sFSb|197+C6p5!Dc49x$BiC,kl!PA@',6#iw9QN%m<\""
        ':5^)/=t8Nol*"tw}I_/;3;)r--.bHW&vcQort\\Q=,SK)VvpSE<7PBC7Bry^aj_4+9_Vr?m1'
        "A5H:OBf.qn4Q|T&fIJoIDLuyc_QNo?d1ZlR2RaD]]dgg8GOpwMLa3we<Idh*Oh2bMDSQz^FS"
        '~\\0|%"(dhu3bB3Em(8^c~X:E.?.oT&F3IORQ@:@1P8\\{55VPkmMgM{nL~81(z^&Eh|dJt0{7'
        "wlwyx#(:Rk=X;qOi[H}QY$>CIOo|bLTO,ZNvsAgQ)Ed^@5\\']E/[z1\\sj/>Xp\\!w+~I1}HV"
        "IG]`,QNts':9M5$E(9sipqGt=,MPgFNO[@iS\"90Y<)c@yw5FSt'E=7@^Z]b+ehEDxO<]ks%cAO"
        "`i{+/(v5.hkku1?EzU1e)AaMvAZLy8ySiz-BoijE}rS<5%A*YkG,wH5g2Qd@*N]Ov#J5ZP%Q$76"
        'G{NMe@|/SKU^8b_Gu}dvj1X56.Mt^u%W~rSwuL.UI-&wX!?/&9LCF4nJS<FRHL"lm<GQ+7a^He'
        "ecsC:Fx.uqXhFum$knu>Mk0OP/'p%)h(nxy?/i;/t!N\"g62Wj(akE5!n719B<&SKqXbt|7O6fZ"
        "~kns3A6H:9k^&3gb1SPO<5Amh348sRc@g#?7[kBJ.%-(,W}P648Obl"
    )
    assert num_alignments == 4100
    assert first_alignment == Result(
        channel=35,
        read_number=3820,
        read_id="0004378f-19e1-4c9e-af1f-850535a43256",
        seq=s,
        decision=Decision.no_seq,
        barcode=None,
        basecall_data=FastqRecord(
            name="0004378f-19e1-4c9e-af1f-850535a43256",
            description="runid=710c3a8095a9c2ea5c598018a4315cebd038270d read=3820 ch=35"
            " start_time=2023-08-16T10:09:15+00:00 flow_cell_id=FAS47120"
            " sample_id=HEK293_LSK110_nemo",
            sequence=s,
            quality=q,
            comment="+",
        ),
        alignment_data=None,
    )


def test_get_fq(tmpdir):
    # Create temporary directory for testing
    directory = tmpdir.mkdir("test_dir")

    # Create some sample fastq files with different extensions
    directory.join("sample1.fastq").write("Sample content")
    directory.join("sample2.fastq.gz").write("Sample content")
    directory.join("sample3.fq").write("Sample content")
    directory.join("sample4.fq.gz").write("Sample content")
    directory.join("sample5.txt").write("Not a fastq file")
    directory.join("sample6.txt.gz").write("Ahhhhhhhhh")
    directory.join("sample7.fastz.xz").write("Wrong compression")

    # Get the list of fastq files using the function being tested
    fastq_files = list(get_fq(directory))

    # Assert that the function returns the correct number of files
    assert len(fastq_files) == 4

    # Assert that all the expected fastq files are in the list
    assert str(directory.join("sample1.fastq")) in fastq_files
    assert str(directory.join("sample2.fastq.gz")) in fastq_files
    assert str(directory.join("sample3.fq")) in fastq_files
    assert str(directory.join("sample4.fq.gz")) in fastq_files

    # Assert that a non-fastq file is not included in the list
    assert str(directory.join("sample5.txt")) not in fastq_files
