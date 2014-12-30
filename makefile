BCALM_BIN_DIR = bin/bcalm
BWA_BIN_DIR=bin/bwa
DALIGNER_BIN_DIR=bin/daligner
HMMEST_BIN_DIR = bin/hmmest
KMC_BIN_DIR = bin/kmc
MINIA_BIN_DIR=bin/minia
SAMTOOLS_BIN_DIR=bin/samtools
RAPSEARCH_BIN_DIR=bin/rapsearch
SAMTOOLS_BIN_DIR=bin/samtools
SOAP_BIN_DIR=bin/soapdenovo
SPADES_BIN_DIR = bin/spades
TRINITY_BIN_DIR = bin/trinity


all: 
	-mkdir -p $(KMC_BIN_DIR) $(BCALM_BIN_DIR) $(HMMEST_BIN_DIR) $(SAMTOOLS_BIN_DIR) $(SOAP_BIN_DIR) $(RAPSEARCH_BIN_DIR) $(DALIGNER_BIN_DIR) $(BWA_BIN_DIR) $(MINIA_BIN_DIR) $(SPADES_BIN_DIR) $(TRINITY_BIN_DIR)
	cd src/bcalm && make && cp bcalm ../../$(BCALM_BIN_DIR)
	cd src/bwa && make && cp bwa ../../$(BWA_BIN_DIR)
	cd src/DALIGNER/DALIGNER && make
	cd src/DALIGNER/DAZZ_DB && make
	cd src/KMC && make && cp bin/* ../../$(KMC_BIN_DIR)
	cd src/HMMEST && make 
	cd src/Minia && make
	cd src/Trinity && make

clean:
	cd src/bcalm && make clean
	cd src/bwa && make clean
	-rm -rf bin
