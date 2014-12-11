#!/bin/bash
#simple assembly test of a synthetic 10 K genome, to verify not completely broken

rm -f t10.contigs.fa

./minia test/read50x_ref10K_e001.fasta 25 3 10000 t10 &> /dev/null


diff t10.contigs.fa  test/result10K.fasta > /dev/null

var=$?

if [ $var -eq 0 ] 
then
    echo Test PASSED
    exit 0
else
    echo Test FAILED
    exit 1
fi
