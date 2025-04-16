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
