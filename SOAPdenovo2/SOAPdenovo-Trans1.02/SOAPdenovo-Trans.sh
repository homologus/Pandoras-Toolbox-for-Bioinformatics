SOAPdenovo-Trans pregraph -s config -K 31 -o output_prefix  
SOAPdenovo-Trans contig -g output_prefix  -e 2 -M 2
SOAPdenovo-Trans map -s config -g output_prefix -r
SOAPdenovo-Trans scaff -g output_prefix -F -L 100 -r  
