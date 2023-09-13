import click


@click.group()
@click.option("--debug/--no-debug")
def cli():
    pass


if __name__ == "__main__":
    cli()
