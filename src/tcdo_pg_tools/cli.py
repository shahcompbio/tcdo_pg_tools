import click
from tcdo_pg_tools.fusion_merge import fusion_merge
from tcdo_pg_tools.coverage_calculator import coverage_calculator
from tcdo_pg_tools.merge_proteome import merge_pg_results, merge_fasta
import importlib.metadata

@click.group()
@click.version_option(version=importlib.metadata.version('tcdo-pg-tools'),
                      message="%(version)s")
def cli():
    pass


cli.add_command(fusion_merge)
cli.add_command(coverage_calculator)
cli.add_command(merge_pg_results)
cli.add_command(merge_fasta)

if __name__ == "__main__":
    cli()
