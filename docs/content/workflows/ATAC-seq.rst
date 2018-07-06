.. _ATAC-seq:

ATAC-seq
============

Input requirements
---------------------------

The DNA mapping pipeline generates output that is fully compatible with the ChIP-seq pipeline input requirements!
When running the ChIP-seq pipeline, please specify the output directory of DNA-mapping pipeline as the working directory (-w).


.. argparse::
    :ref: snakePipes.ATAC-seq.parse_args
    :prog: ATAC-seq
