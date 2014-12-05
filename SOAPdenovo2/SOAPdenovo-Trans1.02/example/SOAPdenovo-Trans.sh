../SOAPdenovo-Trans-31kmer pregraph -s config -K 31 -o output >log  
../SOAPdenovo-Trans-31kmer contig -g output  -e 3 -M 3 >>log
../SOAPdenovo-Trans-31kmer map -s config -g output -r >>log
../SOAPdenovo-Trans-31kmer scaff -g output -F -L 50 -r >>log
