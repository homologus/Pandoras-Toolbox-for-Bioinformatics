KMC_BIN_DIR = bin/kmc
BLAST_BIN_DIR = bin/blast
HMMEST_BIN_DIR = bin/hmmest
SAMTOOLS_BIN_DIR=bin/samtools
SOAP_BIN_DIR=bin/soapdenovo
RAPSEARCH_BIN_DIR=bin/rapsearch
DALIGNER_BIN_DIR=bin/daligner
BWA_BIN_DIR=bin/bwa
MINIA_BIN_DIR=bin/minia
SAILFISH_BIN_DIR=bin/sailfish
SPADES_BIN_DIR = bin/spades
TRINITY_BIN_DIR = bin/trinity


all: 
	-mkdir -p $(KMC_BIN_DIR) $(BLAST_BIN_DIR) $(HMMEST_BIN_DIR) $(SAMTOOLS_BIN_DIR) $(SOAP_BIN_DIR) $(RAPSEARCH_BIN_DIR) $(DALIGNER_BIN_DIR) $(BWA_BIN_DIR) $(MINIA_BIN_DIR) $(SAILFISH_BIN_DIR) $(SPADES_BIN_DIR) $(TRINITY_BIN_DIR)

clean:
	-rm -rf bin

