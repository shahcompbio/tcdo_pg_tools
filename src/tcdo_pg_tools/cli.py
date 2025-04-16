import click
from tcdo_pg_tools.fusion_merge import fusion_merge
from tcdo_pg_tools.coverage_calculator import coverage_calculator

@click.group()
def cli():
    pass


cli.add_command(fusion_merge)
cli.add_command(coverage_calculator)

if __name__ == "__main__":
    cli()
