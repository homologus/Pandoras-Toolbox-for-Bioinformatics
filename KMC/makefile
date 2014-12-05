all: kmc
	
BOOST_LIB = /boost/boost_1_55_0/stage/lib
BOOST_H = /boost/boost_1_55_0

KMC_BIN_DIR = bin
KMC_MAIN_DIR = kmer_counter
KMC_API_DIR = kmc_api
KMC_DUMP_DIR = kmc_dump

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -static -fopenmp -std=c++11 -I $(BOOST_H)
CLINK	= -lm -static -fopenmp -O3 -std=c++11 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

kmc: $(KMC_MAIN_DIR)/kmer_counter.o $(KMC_MAIN_DIR)/mmer.o $(KMC_MAIN_DIR)/mem_disk_file.o  $(KMC_MAIN_DIR)/rev_byte.o $(KMC_MAIN_DIR)/fastq_reader.o $(KMC_MAIN_DIR)/timer.o $(KMC_MAIN_DIR)/radix.o $(KMC_MAIN_DIR)/kb_completer.o $(KMC_MAIN_DIR)/kb_storer.o $(KMC_MAIN_DIR)/kmer.o
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $(KMC_MAIN_DIR)/kmer_counter.o $(KMC_MAIN_DIR)/mem_disk_file.o $(KMC_MAIN_DIR)/rev_byte.o $(KMC_MAIN_DIR)/mmer.o $(KMC_MAIN_DIR)/fastq_reader.o $(KMC_MAIN_DIR)/timer.o $(KMC_MAIN_DIR)/radix.o $(KMC_MAIN_DIR)/kb_completer.o $(KMC_MAIN_DIR)/kb_storer.o $(KMC_MAIN_DIR)/kmer.o $(KMC_MAIN_DIR)/libs/alibelf64.a $(KMC_MAIN_DIR)/libs/libz.a $(KMC_MAIN_DIR)/libs/libbz2.a $(BOOST_LIB)/libboost_thread.a $(BOOST_LIB)/libboost_filesystem.a $(BOOST_LIB)/libboost_system.a

kmc_dump: $(KMC_DUMP_DIR)/nc_utils.o $(KMC_API_DIR)/mmer.o $(KMC_DUMP_DIR)/kmc_dump.o $(KMC_API_DIR)/kmc_file.o $(KMC_API_DIR)/kmer_api.o
	-mkdir -p $(KMC_BIN_DIR)
	$(CC) $(CLINK) -o $(KMC_BIN_DIR)/$@ $(KMC_DUMP_DIR)/nc_utils.o $(KMC_API_DIR)/mmer.o $(KMC_DUMP_DIR)/kmc_dump.o $(KMC_API_DIR)/kmc_file.o $(KMC_API_DIR)/kmer_api.o

clean:
	-rm $(KMC_MAIN_DIR)/*.o
	-rm $(KMC_API_DIR)/*.o
	-rm $(KMC_DUMP_DIR)/*.o
	-rm -rf bin

all: kmc kmc_dump
