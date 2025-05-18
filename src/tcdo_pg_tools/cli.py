import click
from tcdo_pg_tools.fusion_merge import fusion_merge
from tcdo_pg_tools.coverage_calculator import coverage_calculator
from tcdo_pg_tools.fasta_merge import fasta_merge

@click.group()
def cli():
    pass


cli.add_command(fusion_merge)
cli.add_command(coverage_calculator)
cli.add_command(fasta_merge)

if __name__ == "__main__":
    cli()
