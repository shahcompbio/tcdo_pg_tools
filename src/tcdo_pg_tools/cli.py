import click
from tcdo_pg_tools.fusion_merge import fusion_merge

@click.group()
def cli():
    pass


cli.add_command(fusion_merge)

if __name__ == "__main__":
    cli()
