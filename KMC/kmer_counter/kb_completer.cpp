#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.0
  Date   : 2014-07-04
*/
#include <algorithm>
#include <numeric>
#include <iostream>
#include "kb_completer.h"

using namespace std;

extern uint64 total_reads;



//************************************************************************************************************
// CKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Assign queues and monitors
CKmerBinCompleter::CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues)
{
	mm		       = Queues.mm;
	file_name      = Params.output_file_name;
	kq             = Queues.kq;
	bd		       = Queues.bd;
	s_mapper	   = Queues.s_mapper;
	memory_bins    = Queues.memory_bins;

	kmer_file_name = file_name + ".kmc_suf";
	lut_file_name  = file_name + ".kmc_pre";

	kmer_len       = Params.kmer_len;
	signature_len  = Params.signature_len;

	cutoff_min     = Params.cutoff_min;
	cutoff_max     = Params.cutoff_max;
	counter_max    = Params.counter_max;
	lut_prefix_len = Params.lut_prefix_len;

	kmer_t_size    = Params.KMER_T_size;

	use_quake      = Params.use_quake;
}

//----------------------------------------------------------------------------------
CKmerBinCompleter::~CKmerBinCompleter()
{
}

//----------------------------------------------------------------------------------
// Store sorted and compacted bins to the output file
void CKmerBinCompleter::ProcessBins()
{
	int32 bin_id;
	uchar *data = NULL;
	uint64 data_size = 0;
	uchar *lut = NULL;
	uint64 lut_size = 0;
	uint64 counter_size = 0;
	uint32 sig_map_size = (1 << (signature_len * 2)) + 1;
	uint32 *sig_map = new uint32[sig_map_size];
	fill_n(sig_map, sig_map_size, 0);
	uint32 lut_pos = 0;
	if(use_quake)
		counter_size = 4;
	else
		counter_size = min(BYTE_LOG(cutoff_max), BYTE_LOG(counter_max));	
	
	// Open output file
	FILE *out_kmer = fopen(kmer_file_name.c_str(), "wb");
	if(!out_kmer)
	{
		cout << "Error: Cannot create " << kmer_file_name << "\n";
		exit(1);
		return;
	}

	FILE *out_lut = fopen(lut_file_name.c_str(), "wb");
	if(!out_lut)
	{
		cout << "Error: Cannot create " << lut_file_name << "\n";
		fclose(out_kmer);
		exit(1);
		return;
	}

	uint64 _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total;
	uint64 n_recs = 0;

	_n_unique = _n_cutoff_min = _n_cutoff_max = _n_total = 0;
	n_unique  = n_cutoff_min  = n_cutoff_max  = n_total  = 0;

	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";

	// Markers at the beginning
	fwrite(s_kmc_pre, 1, 4, out_lut);
	fwrite(s_kmc_suf, 1, 4, out_kmer);

	// Process priority queue of ready-to-output bins
	while(!kq->empty())
	{
		// Get the next bin
		if (!kq->pop(bin_id, data, data_size, lut, lut_size, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total))
			continue;

		// Decrease memory size allocated by stored bin
		string name;
		uint64 n_rec;
		uint64 n_plus_x_recs;
		uint64 n_super_kmers;
		uint64 raw_size;
		CMemDiskFile *file;

		bd->read(bin_id, file, name, raw_size, n_rec, n_plus_x_recs, n_super_kmers);

		uint64 lut_recs        = lut_size / sizeof(uint64);
		
		// Write bin data to the output file
		fwrite(data, 1, data_size, out_kmer);
		memory_bins->free(bin_id, CMemoryBins::mba_suffix);

		uint64 *ulut = (uint64*) lut;
		for(uint64 i = 0; i < lut_recs; ++i)
		{
			uint64 x  = ulut[i];
			ulut[i]   = n_recs;
			n_recs   += x;
		}
		fwrite(lut, lut_recs, sizeof(uint64), out_lut);
		//fwrite(&n_rec, 1, sizeof(uint64), out_lut);
		memory_bins->free(bin_id, CMemoryBins::mba_lut);

		n_unique	 += _n_unique;
		n_cutoff_min += _n_cutoff_min;
		n_cutoff_max += _n_cutoff_max;
		n_total      += _n_total;
		for (uint32 i = 0; i < sig_map_size; ++i)
		{
			if (s_mapper->get_bin_id(i) == bin_id)
			{
				sig_map[i] = lut_pos;
			}
		}
		++lut_pos;
	}
	
	// Marker at the end
	fwrite(s_kmc_suf, 1, 4, out_kmer);
	fclose(out_kmer);

	fwrite(&n_recs, 1, sizeof(uint64), out_lut);

	//store signature mapping 
	fwrite(sig_map, sizeof(uint32), sig_map_size, out_lut);	

	// Store header
	uint32 offset = 0;

	store_uint(out_lut, kmer_len, 4);				offset += 4;
	store_uint(out_lut, (uint32) use_quake, 4);		offset += 4;	// mode: 0 (counting), 1 (Quake-compatibile counting)
	store_uint(out_lut, counter_size, 4);			offset += 4;
	store_uint(out_lut, lut_prefix_len, 4);			offset += 4;
	store_uint(out_lut, signature_len, 4);			offset += 4; 
	store_uint(out_lut, cutoff_min, 4);				offset += 4;
	store_uint(out_lut, cutoff_max, 4);				offset += 4;
	store_uint(out_lut, n_unique - n_cutoff_min - n_cutoff_max, 8);		offset += 8;

	// Space for future use
	for(int32 i = 0; i < 7; ++i)
	{
		store_uint(out_lut, 0, 4);
		offset += 4;
	}
	
	store_uint(out_lut, 0x200, 4);
	offset += 4;

	store_uint(out_lut, offset, 4);

	// Marker at the end
	fwrite(s_kmc_pre, 1, 4, out_lut);
	fclose(out_lut);
	cout << "\n";

	delete[] sig_map;
}

//----------------------------------------------------------------------------------
// Return statistics
void CKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	_n_unique	  = n_unique;
	_n_cutoff_min = n_cutoff_min;
	_n_cutoff_max = n_cutoff_max;
	_n_total      = n_total;
}

//----------------------------------------------------------------------------------
// Store single unsigned integer in LSB fashion
bool CKmerBinCompleter::store_uint(FILE *out, uint64 x, uint32 size)
{
	for(uint32 i = 0; i < size; ++i)
		putc((x >> (i * 8)) & 0xFF, out);

	return true;
}


//************************************************************************************************************
// CWKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CWKmerBinCompleter::CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues)
{
	kbc = new CKmerBinCompleter(Params, Queues);
}

//----------------------------------------------------------------------------------
// Destructor
CWKmerBinCompleter::~CWKmerBinCompleter()
{
	delete kbc;
}

//----------------------------------------------------------------------------------
// Execution
void CWKmerBinCompleter::operator()()
{
	kbc->ProcessBins();
}

//----------------------------------------------------------------------------------
// Return statistics
void CWKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	if(kbc)
		kbc->GetTotal(_n_unique, _n_cutoff_min, _n_cutoff_max, _n_total);
}

// ***** EOF
