#!/bin/bash -ve

../../Trinity --seqType fq --single reads.left.fq.gz  --SS_lib_type R  --CPU 2 --JM 1G -o trin_SE_outdir --trimmomatic --normalize_reads



