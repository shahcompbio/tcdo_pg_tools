import click
from tcdo_pg_tools.merge_pg_results import merge_pg_results

@click.command()
@click.option('-i', '--input_csv', required=True, type=click.Path(exists=True),
              help='three column csv (fasta: fasta path, '
                   'name: sample name, condition: condition)')
@click.option('-t', '--info_table', required=False,
              default='info_table.tsv',
              type=click.Path(), help="Path to index tsv for merged protein IDs")
@click.option('-fa','--merged_fasta', required=False,
              type=click.Path(), default='merged.fasta',
              help="Path to merged fasta file")
@click.option('--upset', is_flag=True, default=False, help="plot upset")
@click.option('--upset_path', required=False,
              type=click.Path(), default='upset_plot.svg',
              help="Path to upset plot")
def merge_fasta(input_csv, info_table, merged_fasta, upset, upset_path):
    """
    merge multiple fasta on sequence identity
    """
    merge_pg_results(input_csv, info_table, merged_fasta, upset, upset_path, unique_proteins=False)
    return
