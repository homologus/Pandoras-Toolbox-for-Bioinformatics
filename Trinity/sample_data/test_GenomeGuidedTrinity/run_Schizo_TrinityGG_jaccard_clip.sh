#!/bin/bash -ve

../../Trinity --seqType fq --left top100k.Left.fq.gz --right top100k.Right.fq.gz --genome top100k.genome.gz --genome_guided_use_bam SP2.chr.bam --JM 1G --CPU 2 --genome_guided_max_intro 1000 --jaccard_clip --SS_lib_type RF --output test_Schizo_GG_jaccard_RF_outdir
