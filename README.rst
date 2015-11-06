Identifying high resolution reference points using ChIP-exo
================

Introduction
-------------

The pipeline proposed here can be used to find high resolution reference coordinates.The pipeline also gives the following:
1) Genomic regions containing over represented combination of factors.
2) Over-represented motifs in these enriched regions.
3) Assessment of each candidate reference as a potential organizing center.
4) Distance of reference features from each other.


Requirements
------------
- Python 2.7 or higher.
- Python modules: pysam, pybedtools, numpy, scipt, matplotlib.

Input
-------

- A directory containing peak calls in gff format for all the factors.
- A directory containing indexed BAM files for the same factors.
- Both peak filenames and BAM filenames should start with factorname followed by "_".


Installing and Running the scripts
-----------------------------------
To run the pipeline type the following::
    $ python Run_pipeline.py -h
    $ Usage: Run_pipeline.py [options] input_paths

    $ input_paths may be:
    $ - a directory to run on all files in them

    $ example usages:
    $ python Run_pipeline.py [OPTIONS]
    $ **** Remember that your BAM file and your peak file should start with "factorname_"

Options:
  -h, --help  show this help message and exit
  -p PEAKDIR  The directory containing the peaks call files in gff format.
  -i BAMDIR   The directory containing the BAM files.
  -w WINDOW   Window size, default = 20
  -g GFILE    File containing the chromosome number and length.
  -s GSIZE    Mappable genome size: sg11 = 11,332,237 (default),
              hg19=248,988,565, mm9=2,178,433,024 for read length = 36. Refer
              PMID:22276185
  -v PVAL     P-value cutoff for significant enrichment over background,
              default = 0.05
  -l ILEN     Length of chromHMM segmentation interval, default = 200
  -n STATE    Number of states, default = 12
  -r FASTA    Reference FASTA file.
  -t TSS      File containing TSS coordinates.
  -u UP       Upstream distance from candidate reference point, default = 500
  -d DOWN     Downstream distance from candidate reference point, default =
              500


Help menu
-----------



Output
------



 

.. _Python: https://www.python.org/
.. _pysam: https://code.google.com/p/pysam/
.. _pybedtools: https://pythonhosted.org/pybedtools/
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _matplotlib: http://matplotlib.org/
.. _gff: http://genome.ucsc.edu/FAQ/FAQformat#format3
.. _BAM: https://samtools.github.io/hts-specs/SAMv1.pdf
