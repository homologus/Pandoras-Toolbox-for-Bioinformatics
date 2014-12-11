CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC=Pool.cpp Bank.cpp Bloom.cpp Hash16.cpp LargeInt.cpp Kmer.cpp Terminator.cpp Traversal.cpp LinearCounter.cpp Set.cpp Utils.cpp SortingCount.cpp Debloom.cpp OAHash.cpp
EXEC=minia
OBJ= $(SRC:.cpp=.o)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif

ifeq ($(unitig),1)
 CFLAGS+= -DUNITIG
endif


ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

ifeq ($(omp),1)
 CFLAGS=-O4  -fopenmp -DOMP=1
endif


k := 0$(k) # dummy k if not specified 
K_BELOW_32 := $(shell echo $(k)\<=32 | bc)
K_BELOW_64 := $(shell echo $(k)\<=64 | bc)
ARCH := $(shell getconf LONG_BIT) # detects sizeof(int)
USING_UINT128 := 0
largeintlib := 0

ifeq ($(K_BELOW_32),0)

    # use uint128 when k<=64 and 64-bit architecture
    ifeq ($(K_BELOW_64),1)
        ifeq ($(strip $(ARCH)),64)
            CFLAGS += -Dkmer_type=__uint128_t
            USING_UINT128 := 1
        endif
    endif
    
	# use a bigint library otherwise
    ifeq ($(USING_UINT128),0)
        largeintlib := largeint#ttmath
    endif
endif

# ttmath (now, largeint) is used when you type "make k=[kmer size]" with a kmer size longer than supported integer type,
ifeq ($(largeintlib),ttmath)
    KMER_PRECISION := $(shell echo \($(k)+15\)/16 | bc)
endif
ifeq ($(largeintlib),largeint)
    KMER_PRECISION := $(shell echo \($(k)+31\)/32 | bc)
endif
ifneq ($(largeintlib),0)
    CFLAGS += -D_$(largeintlib) -DKMER_PRECISION=$(KMER_PRECISION)
endif

all:
	$(MAKE) clean
	$(MAKE) $(EXEC)

minia: $(OBJ) Minia.cpp
	$(CXX) -o $@ $(OBJ) Minia.cpp $(CFLAGS) -lz

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)


%.o: %.c %.h 
	$(CXX) -o $@ -c $< $(CFLAGS)
    
clean:
	rm -rf *.o

install:
	cp minia /usr/local/bin
