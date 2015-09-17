Discovering Cis-Regulatory Modules using ChIP-exo
================

Introduction
-------------

The pipeline here can be used to find significant over-represented combination of transcription factors in the genome. The pipeline also gives the following:
1) Information about enriched states (combinations of factors).
2) Over-represented motifs in all enriched states.
3) Distance of reference features from enriched states.
4) Information about different reference points that might be the organization centers of these modules.


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
