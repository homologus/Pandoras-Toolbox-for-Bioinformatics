BCALM_BIN_DIR = bin/bcalm
BWA_BIN_DIR=bin/bwa
DALIGNER_BIN_DIR=bin/daligner
HMMEST_BIN_DIR = bin/hmmest
KMC_BIN_DIR = bin/kmc
MINIA_BIN_DIR=bin/minia
SAMTOOLS_BIN_DIR=bin/samtools
RAPSEARCH_BIN_DIR=bin/rapsearch
SOAP_BIN_DIR=bin/soapdenovo
SPADES_BIN_DIR = bin/spades
TRINITY_BIN_DIR = bin/trinity


all: 
	-mkdir -p $(BCALM_BIN_DIR) $(BWA_BIN_DIR) $(DALIGNER_BIN_DIR) $(HMMEST_BIN_DIR) $(KMC_BIN_DIR) $(MINIA_BIN_DIR) $(SAMTOOLS_BIN_DIR) $(SOAP_BIN_DIR) $(RAPSEARCH_BIN_DIR) $(SPADES_BIN_DIR) $(TRINITY_BIN_DIR)
	cd src/bcalm && make && cp bcalm ../../$(BCALM_BIN_DIR)
	cd src/bwa && make && cp bwa ../../$(BWA_BIN_DIR)
	cd src/DALIGNER/DALIGNER && make
	cd src/DALIGNER/DAZZ_DB && make
	cd src/HMMEST && make 
	cd src/KMC && make && cp bin/* ../../$(KMC_BIN_DIR)
	cd src/Minia && make
	cd src/RAPSearch2 && make
	cd src/samtools/samtools-1.1  && make
	cd src/samtools/bcftools-1.1  && make
	cd src/samtools/htslib-1.1  && make
	cd src/SOAPdenovo2/SOAPdenovo2-src-r240 && make
	cd src/SOAPdenovo2/SOAPdenovo-Trans-src-v1.04 && sh make.sh
	cd src/Trinity && make
	cd src/SPAdes && sh spades_compile.sh

clean:
	cd src/bcalm && make clean
	cd src/bwa && make clean
	cd src/bcalm && make clean
	cd src/bwa && make clean
	cd src/DALIGNER/DALIGNER && make clean
	cd src/DALIGNER/DAZZ_DB && make clean
	cd src/HMMEST && make clean
	cd src/KMC && make clean
	cd src/Minia && make clean
	cd src/RAPSearch2 && make clean
	cd src/samtools/samtools-1.1  && make clean
	cd src/samtools/bcftools-1.1  && make clean
	cd src/samtools/htslib-1.1  && make clean
	cd src/SOAPdenovo2/SOAPdenovo2-src-r240 && make clean
	cd src/SOAPdenovo2/SOAPdenovo-Trans-src-v1.04 && sh clean.sh
	cd src/Trinity && make clean
	rm -rf bin
