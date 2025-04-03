=====
Usage
=====

To use tcdo-pg-tools in a project::

    import tcdo_pg_tools

Merging RNA fusions across samples
------------------------

Currently only for results called by CTAT-LR-Fusion or processed using FusionInspector.
and requires metadata sheet from our internal isabl_utils package. Usage::

    fusion_merge --help
    usage: fusion_merge [-h]
                        input_metadata app_version result_type fusion_table
                        output_fasta

    merge fusion calls across multiple samples based on gene symbols and
    breakpoints.

    positional arguments:
      input_metadata  metadata sheet from isabl_utils w/ paths to FusionInspector
                      results
      app_version     filter for results from specific isabl app version
      result_type     filter on result type in metadata sheet
      fusion_table    path to output fusion info table
      output_fasta    path to fasta output
