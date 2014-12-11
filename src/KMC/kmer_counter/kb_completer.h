/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.0
  Date   : 2014-07-04
*/
#ifndef _KB_COMPLETER_H
#define _KB_COMPLETER_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include "radix.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <stdio.h>


//************************************************************************************************************
// CKmerBinCompleter - complete the sorted bins and store in a file
//************************************************************************************************************
class CKmerBinCompleter {
	CMemoryMonitor *mm;
	string file_name, kmer_file_name, lut_file_name;
	CKmerQueue *kq;
	CBinDesc *bd;
	CSignatureMapper *s_mapper;

	CMemoryBins *memory_bins;
	uint32 lut_prefix_len;
	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total;
	uint32 kmer_t_size;
	int32 cutoff_min, cutoff_max;
	int32 counter_max;
	int32 kmer_len;
	int32 signature_len;
	bool use_quake;

	bool store_uint(FILE *out, uint64 x, uint32 size);

public:
	CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CKmerBinCompleter();

	void ProcessBins();
	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
};


//************************************************************************************************************
// CWKmerBinCompleter - wrapper for multithreading purposes
//************************************************************************************************************
class CWKmerBinCompleter {
	CKmerBinCompleter *kbc;

public:
	CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues);
	~CWKmerBinCompleter();

	void operator()();

	void GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total);
};

#endif

// ***** EOF
