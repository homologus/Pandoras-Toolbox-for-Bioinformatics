BCALM_BIN_DIR = bin/bcalm
BWA_BIN_DIR=bin/bwa
DALIGNER_BIN_DIR=bin/daligner
HMMEST_BIN_DIR = bin/hmmest
KMC_BIN_DIR = bin/kmc
MINIA_BIN_DIR=bin/minia
SAMTOOLS_BIN_DIR=bin/samtools
RAPSEARCH_BIN_DIR=bin/rapsearch
SAILFISH_BIN_DIR=bin/sailfish
SOAP_BIN_DIR=bin/soapdenovo
SPADES_BIN_DIR = bin/spades
TRINITY_BIN_DIR = bin/trinity

DALIGNER_FILES=daligner HPCdaligner LAsort LAmerge LAshow LAsplit LAcat LAcheck
DAZZ_FILES=fasta2DB DB2fasta quiva2DB DB2quiva DBsplit DBdust Catrack DBshow DBstats DBrm simulator
HMM_FILES=alimask hmmalign hmmc2 hmmconvert hmmemit hmmfetch hmmlogo hmmpgmd hmmpgmd_client_example.pl hmmpress hmmpress.itest.pl hmmscan hmmsearch hmmsim hmmstat jackhmmer nhmmer nhmmscan phmmer
CHRYSALIS_FILES=chrysalis.notes QuantifyGraph GraphFromFasta ReadsToComponents.pl IsoformAugment ReadsToTranscripts BreakTransByPairs  JoinTransByPairs  RunButterfly checkLock Chrysalis TranscriptomeFromVaryK


all: 
	-rm -rf bin
	-mkdir -p $(BCALM_BIN_DIR) $(BWA_BIN_DIR) $(DALIGNER_BIN_DIR) $(DALIGNER_BIN_DIR)/Daligner $(DALIGNER_BIN_DIR)/db $(HMMEST_BIN_DIR) $(KMC_BIN_DIR) $(MINIA_BIN_DIR) $(SAMTOOLS_BIN_DIR) $(SOAP_BIN_DIR) $(RAPSEARCH_BIN_DIR) $(SPADES_BIN_DIR) $(TRINITY_BIN_DIR) $(SAILFISH_BIN_DIR)  $(TRINITY_BIN_DIR)/Inchworm $(TRINITY_BIN_DIR)/Chrysalis $(TRINITY_BIN_DIR)/Butterfly

	cd src/bcalm && make && cp bcalm ../../$(BCALM_BIN_DIR)
	cd src/bwa && make && cp bwa ../../$(BWA_BIN_DIR)
	cd src/DALIGNER/DALIGNER && make  && cp $(DALIGNER_FILES) ../../../$(DALIGNER_BIN_DIR)/Daligner 
	cd src/DALIGNER/DAZZ_DB && make && cp $(DAZZ_FILES) ../../../$(DALIGNER_BIN_DIR)/db
	cd src/HMMEST && ./configure && make && cd src && cp $(HMM_FILES) ../../../$(HMMEST_BIN_DIR) 
	cd src/KMC && make && cp bin/* ../../$(KMC_BIN_DIR)
	cd src/Minia && make  && cp minia ../../$(MINIA_BIN_DIR)
	cd src/RAPSearch2 && rm -f boost lib*a && make clean && ln -s ../../boost_1_55_0/boost boost && ln -s ../../boost_1_55_0/stage/lib/libboost_chrono.a libboost_chrono.a && ln -s ../../boost_1_55_0/stage/lib/libboost_serialization.a libboost_serialization.a && ln -s ../../boost_1_55_0/stage/lib/libboost_system.a libboost_system.a && ln -s ../../boost_1_55_0/stage/lib/libboost_thread.a libboost_thread.a && make && cp *search ../../$(RAPSEARCH_BIN_DIR)
	cd src/samtools/samtools-1.1  && make && cp samtools ../../../$(SAMTOOLS_BIN_DIR) 
	cd src/samtools/bcftools-1.1  && make && cp bcftools ../../../$(SAMTOOLS_BIN_DIR) 
	cd src/samtools/htslib-1.1  && make && cp bgzip ../../../$(SAMTOOLS_BIN_DIR)
	cd src/SOAPdenovo2/SOAPdenovo2-src-r240 && make && cp SOAP* ../../../$(SOAP_BIN_DIR)
	cd src/SOAPdenovo2/SOAPdenovo-Trans-src-v1.04 && sh make.sh && cp SOAP* ../../../$(SOAP_BIN_DIR)
	cd src/SPAdes && rm -f ext/include/boost && ln -s ../../../../boost_1_55_0/boost ext/include/boost && sh spades_compile.sh && cp build_spades/bin/* ../../$(SPADES_BIN_DIR) && cp *py ../../$(SPADES_BIN_DIR)
	cd src/Trinity && make && cp Trinity ../../$(TRINITY_BIN_DIR) && cp Inchworm/bin/* ../../$(TRINITY_BIN_DIR)/Inchworm && cd Chrysalis && cp $(CHRYSALIS_FILES) ../../../$(TRINITY_BIN_DIR)/Chrysalis && cd .. && cp Butterfly/Butterfly.jar ../../$(TRINITY_BIN_DIR)/Butterfly

clean:
	cd src/bcalm && make clean
	cd src/bwa && make clean
	cd src/bcalm && make clean
	cd src/bwa && make clean
	cd src/DALIGNER/DALIGNER && make clean
	cd src/DALIGNER/DAZZ_DB && make clean
	cd src/KMC && make clean
	cd src/Minia && make clean
	cd src/RAPSearch2 && make clean 
	cd src/samtools/samtools-1.1  && make clean
	cd src/samtools/bcftools-1.1  && make clean
	cd src/samtools/htslib-1.1  && make clean
	cd src/SOAPdenovo2/SOAPdenovo2-src-r240 && make clean
	cd src/SOAPdenovo2/SOAPdenovo-Trans-src-v1.04 && sh clean.sh
	cd src/SPAdes && rm -rf build_spades
	cd src/Trinity && make clean
