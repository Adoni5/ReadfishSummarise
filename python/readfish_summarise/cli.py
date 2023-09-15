import click
from readfish_summarise.summarise import fastq


@click.group()
# @click.option("--debug/--no-debug")
def cli():
    pass


cli.add_command(fastq)

if __name__ == "__main__":
    cli()
