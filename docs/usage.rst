=====
Usage
=====

To use tcdo-pg-tools in a project::

    import tcdo_pg_tools

See Commandline Interface section for commandline options.


Examples
------
Computing AA coverage across multi-enzyme digests
----
To compute AA seq coverage across different combinations of enzymes, use the `coverage-calculator` feature. For example::

    tcdo_pg_tools coverage-calculator --fragpipe_dir /Volumes/kentsis/proteomics/fragpipe_results/APS010.1_A673_proteogenomics/spike_in --enzymes argc,aspn,gluc,in-house_chymotrypsin,lysc,lysn,proalanase,trypsin --output_tsv protein_coverage.tsv

The fragpipe output directory must contain subdirectories with the result of each enzyme,
and these directories must be labeled starting with the enzyme name (e.g., `trypsin_diaPASEF`)

Merging multiple proteomegenerator fasta
------
To merge across multiple proteomegenerator fasta file by merging proteins
with the same amino acid sequence, use the `merge-fasta` feature::

    tcdo_pg_tools merge-fasta -i input.csv

The input.csv must have three columns: fasta, sample, condition.
`fasta` is the path to the protein fasta file, `sample` is the sample name, and `condition` is the condition for a given sample (e.g., tumor, normal).

You can use the `--upset` flag to output an upset plot likeso::

    tcdo_pg_tools merge-fasta -i input.csv --upset

The upset plot will be plot across the `condition` column.

Merging proteomegenerator results across multiple samples
--------
To merge multiple proteomegenerator results, use `merge-pg-results`, which behaves in the same way
as `merge-fasta`, except it filters on the protein.tsv output of philosopher to identify unique proteins.
So you can run::

    tcdo_pg_tools merge-pg-results -i input.csv --upset

Where here the input.csv must have four columns: fasta, protein_table, sample, condition.
`fasta` is the path to the protein fasta file, `protein_table` is the `protein.tsv` file that is output by Philosopher (in the Fragpipe output directory),
`sample` is the sample name, and `condition` is the condition for a given sample (e.g., tumor, normal).
