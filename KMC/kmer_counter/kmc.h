/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.0
  Date   : 2014-07-04
*/

#ifndef _KMC_H
#define _KMC_H

#include "defs.h"
#include "params.h"
#include "kmer.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <numeric>
#include "queues.h"
#include "timer.h"
#include "fastq_reader.h"
#include "kb_collector.h"
#include "kb_completer.h"
#include "kb_reader.h"
#include "kb_sorter.h"
#include "kb_storer.h"
#include "s_mapper.h"
#include "splitter.h"
#include "libs/asmlib.h"
#include <boost/filesystem.hpp>

#ifdef DEVELOP_MODE
#include "develop.h"
#endif

using namespace std;


template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> class CKMC {
	bool initialized;

	CStopWatch w0, heuristic_time , w1, w2;

	// Parameters (input and internal)
	CKMCParams Params;

	// Memory monitor and queues
	CKMCQueues Queues;

	// Thread groups
	vector<thread> gr0_1, gr0_2;
	vector<thread> gr1_1, gr1_2, gr1_3, gr1_4, gr1_5;		// thread groups for 1st stage
	vector<thread> gr2_1, gr2_2, gr2_3;						// thread groups for 2nd stage

	uint64 n_unique, n_cutoff_min, n_cutoff_max, n_total, n_reads, tmp_size, n_total_super_kmers;

	// Threads
	vector<CWStatsFastqReader*> w_stats_fastqs;
	vector<CWStatsSplitter<false>*> w_stats_splitters;
	vector<CWFastqReader*> w_fastqs;
	vector<CWSplitter<QUAKE_MODE>*> w_splitters;
	CWKmerBinStorer *w_storer;

	CWKmerBinReader<KMER_T, SIZE>* w_reader;
	vector<CWKmerBinSorter<KMER_T, SIZE>*> w_sorters;
	CWKmerBinCompleter *w_completer;

	void SetThreads1Stage();
	void SetThreads2Stage(vector<int64>& sorted_sizes);
	
	bool AdjustMemoryLimits();
	void AdjustMemoryLimitsStage2();

	void ShowSettingsStage1();
	void ShowSettingsStage2();

public:
	CKMC();
	~CKMC();
	
	void SetParams(CKMCParams &_Params);
	bool Process();
	void GetStats(double &time1, double &time2, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64& _n_total_super_kmers);
};


//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> CKMC<KMER_T, SIZE, QUAKE_MODE>::CKMC()
{
// OpenMP support is a must, so do not compile if it is not supported
#if !defined(_OPENMP)
	BOOST_STATIC_ASSERT_MSG(false, "You need to use OpenMP");
#endif

	initialized   = false;
	Params.kmer_len      = 0;
	Params.n_readers     = 1;
	Params.n_splitters   = 1;
	Params.n_sorters     = 1;
	//Params.n_omp_threads = 1;
	Queues.s_mapper = NULL;
}

//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> CKMC<KMER_T, SIZE, QUAKE_MODE>::~CKMC()
{
}

//----------------------------------------------------------------------------------
// Set params of the k-mer counter
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::SetParams(CKMCParams &_Params)
{
	Params = _Params;
	Params.kmer_len	= Params.p_k;

	if (Params.kmer_len % 32 == 0)
		Params.max_x = 0;
	else
		Params.max_x = MIN(31 - (Params.kmer_len % 32), KMER_X);

	Params.verbose	= Params.p_verbose;	
	// Technical parameters related to temporary files
	
	Params.signature_len	 = Params.p_p1;
	Params.bin_part_size     = 1 << 16; 
	
	
	// Thresholds for counters
	Params.cutoff_min   = Params.p_ci;
	Params.cutoff_max   = Params.p_cx;
	Params.counter_max  = Params.p_cs;
	Params.use_quake    = Params.p_quake;

	Params.lowest_quality = Params.p_quality;
	Params.both_strands   = Params.p_both_strands;
	Params.mem_mode		  = Params.p_mem_mode;
	
	// Technical parameters related to no. of threads and memory usage
	if(Params.p_sf && Params.p_sp && Params.p_so && Params.p_sr)
	{
		Params.n_readers     = NORM(Params.p_sf, 1, 32);
		Params.n_splitters   = NORM(Params.p_sp, 1, 32);
		Params.n_sorters     = NORM(Params.p_sr, 1, 32);
		//Params.n_omp_threads = NORM(Params.p_so, 1, 32);
		Params.n_omp_threads.assign(Params.n_sorters, NORM(Params.p_so, 1, 32));
	}
	else
	{
		// Adjust the number of threads according to the current hardware
		Params.n_threads = Params.p_t;
		if (!Params.n_threads)
			Params.n_threads = thread::hardware_concurrency();
		SetThreads1Stage();
	}

	//Params.max_mem_size  = NORM(((uint64) Params.p_m) << 30, (uint64) MIN_MEM << 30, 1024ull << 30);
	Params.max_mem_size = NORM(((uint64)Params.p_m) * 1000000000ull, (uint64)MIN_MEM * 1000000000ull, 1024ull * 1000000000ull);

	Params.file_type		= Params.p_file_type;

	Params.KMER_T_size = sizeof(KMER_T);

	initialized = true; 

	SetMemcpyCacheLimit(8);			// Sets the asmlib's memcpy function to make copy without use of cache memory
}

//----------------------------------------------------------------------------------
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::SetThreads1Stage()
{
	if (!Params.p_sf || !Params.p_sp || !Params.p_sr || !Params.p_so)
	{
		int cores = Params.n_threads;
		bool gz_bz2 = false;
		vector<uint64> file_sizes;
		
		for (auto& p : Params.input_file_names)
		{
			string ext(p.end() - 3, p.end());
			if (ext == ".gz" || ext == ".bz2")
			{
				gz_bz2 = true;
				//break;
			}
			FILE* tmp = my_fopen(p.c_str(), "rb");
			if (!tmp)
			{
				cout << "Cannot open file: " << p.c_str();
				exit(1);
			}
			my_fseek(tmp, 0, SEEK_END);
			file_sizes.push_back(my_ftell(tmp));
			fclose(tmp);
		}
		if (gz_bz2)
		{
			sort(file_sizes.begin(), file_sizes.end(), greater<uint64>());
			uint64 file_size_threshold = (uint64)(file_sizes.front() * 0.05);
			int32 n_allowed_files = 0;
			for(auto& p : file_sizes)
			if (p > file_size_threshold)
				++n_allowed_files;
			Params.n_readers = MIN(n_allowed_files, MAX(1, cores / 2));
		}
		else
			Params.n_readers = 1;
		Params.n_splitters = MAX(1, cores - Params.n_readers);
	}
}
//----------------------------------------------------------------------------------
template<typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::SetThreads2Stage(vector<int64>& sorted_sizes)
{	
	if (!Params.p_sf || !Params.p_sp || !Params.p_sr || !Params.p_so)
	{
		if (Params.n_threads == 1)
		{
			Params.n_sorters = 1;
			Params.n_omp_threads.assign(1, 1);
		}
		else
		{			
			int64 _10th_proc_bin_size = MAX(sorted_sizes[int(sorted_sizes.size() * 0.1)], 1);			
			Params.n_sorters = NORM(Params.max_mem_size / _10th_proc_bin_size, 1, Params.n_threads);
			Params.n_omp_threads.assign(Params.n_sorters, MAX(1, Params.n_threads / Params.n_sorters));
			int threads_left = Params.n_threads - Params.n_omp_threads.front() * Params.n_sorters;
			for (uint32 i = 0; threads_left; --threads_left, ++i)
				Params.n_omp_threads[i%Params.n_sorters]++;
		}
	}
}

template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::AdjustMemoryLimitsStage2()
{
	// Memory for 2nd stage
	// Settings for memory manager of radix internal buffers
	Params.mem_part_pmm_radix_buf = (256 * BUFFER_WIDTH + ALIGNMENT) * sizeof(uint64);


	int64 sum_n_omp_threads = 0;
	for (auto& p : Params.n_omp_threads)
		sum_n_omp_threads += p;

	//Params.mem_tot_pmm_radix_buf = Params.mem_part_pmm_radix_buf * Params.n_sorters * Params.n_omp_threads;
	
	Params.mem_tot_pmm_radix_buf = Params.mem_part_pmm_radix_buf * sum_n_omp_threads;


	if (Params.use_quake)
	{
		Params.mem_part_pmm_prob = (CKmerBinSorter<KMER_T, SIZE>::PROB_BUF_SIZE + 1) * sizeof(double);
		Params.mem_tot_pmm_prob = Params.n_sorters * Params.mem_part_pmm_prob;
	}
	else
		Params.mem_part_pmm_prob = Params.mem_tot_pmm_prob = 0;
	if (!Params.use_quake && Params.both_strands)
	{
		Params.mem_part_pmm_epxand = EXPAND_BUFFER_RECS * sizeof(KMER_T);
		Params.mem_tot_pmm_epxand = sum_n_omp_threads * Params.mem_part_pmm_epxand;
	}
	else
		Params.mem_part_pmm_epxand = Params.mem_tot_pmm_epxand = 0;

	Params.max_mem_stage2 = Params.max_mem_size - Params.mem_tot_pmm_radix_buf - Params.mem_tot_pmm_prob - Params.mem_tot_pmm_epxand;
}

//----------------------------------------------------------------------------------
// Adjust the memory limits for queues and other large data structures
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> bool CKMC<KMER_T, SIZE, QUAKE_MODE>::AdjustMemoryLimits()
{
	// Memory for splitter internal buffers
	int64 m_rest = Params.max_mem_size;  

	Params.mem_part_pmm_stats = ((1 << Params.signature_len * 2) + 1) * sizeof(uint32);
	Params.mem_tot_pmm_stats = (Params.n_splitters + 1 + 1) * Params.mem_part_pmm_stats; //1 merged in main thread, 1 for sorting indices

	
	// Settings for memory manager of FASTQ buffers
	Params.fastq_buffer_size = 32 << 20;
	do {
		if(Params.fastq_buffer_size & (Params.fastq_buffer_size-1))
			Params.fastq_buffer_size &= Params.fastq_buffer_size - 1;
		else
			Params.fastq_buffer_size = Params.fastq_buffer_size / 2 + Params.fastq_buffer_size / 4;
		Params.mem_part_pmm_fastq = Params.fastq_buffer_size + CFastqReader::OVERHEAD_SIZE;
		Params.mem_tot_pmm_fastq  = Params.mem_part_pmm_fastq * (Params.n_readers + Params.n_splitters + 96);
	} while(Params.mem_tot_pmm_fastq > m_rest * 0.17);
	m_rest -= Params.mem_tot_pmm_fastq;

	// Subtract memory for buffers for decompression of FASTQ files
	while(Params.n_readers * Params.gzip_buffer_size > m_rest / 10)
		Params.gzip_buffer_size /= 2;
	m_rest -= Params.n_readers * Params.gzip_buffer_size;

	// Subtract memory for bin collectors internal buffers
	m_rest -= Params.n_splitters * Params.bin_part_size * sizeof(KMER_T);

	// Settings for memory manager of reads
	Params.mem_part_pmm_reads = (CSplitter<QUAKE_MODE>::MAX_LINE_SIZE + 1) * sizeof(double);
	Params.mem_tot_pmm_reads  = Params.mem_part_pmm_reads * 2 * Params.n_splitters;
	m_rest -= Params.mem_tot_pmm_reads;

	// Max. memory for single package
	Params.max_mem_storer_pkg = 1ll << 25; 

	Params.mem_part_pmm_bins = Params.bin_part_size;

	Params.mem_tot_pmm_bins = m_rest;

	// memory for storer internal buffer
	if(Params.max_mem_size >= 16ll << 30)
		Params.max_mem_storer = (int64) (Params.mem_tot_pmm_bins * 0.75);
	else
		Params.max_mem_storer = (int64) (Params.mem_tot_pmm_bins * 0.65);

	if(Params.max_mem_storer < (1ll << 28))
		return false;

	return true;
}

//----------------------------------------------------------------------------------
// Show the settings of the KMC (in verbose mode only)
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::ShowSettingsStage1()
{
	if(!Params.verbose)
		return;

	cout << "\n********** Used parameters: **********\n";

	cout << "No. of input files           : " << Params.input_file_names.size() << "\n";
	cout << "Output file name             : " << Params.output_file_name << "\n";
	cout << "No. of working directories   : " << 1 << "\n";
	cout << "Input format                 : "; 
	switch (Params.file_type)
	{
	case fasta:
		cout << "FASTA\n";
		break;
	case fastq:
		cout << "FASTQ\n";
		break;
	case multiline_fasta:
		cout << "MULTI LINE FASTA\n";
		break;
	}
	cout << "\n";
	cout << "k-mer length                 : " << Params.kmer_len << "\n";
	cout << "Max. k-mer length            : " << MAX_K << "\n";
	cout << "Signature length             : " << Params.signature_len << "\n"; 
	cout << "Min. count threshold         : " << Params.cutoff_min << "\n";
	cout << "Max. count threshold         : " << Params.cutoff_max << "\n";
	cout << "Max. counter value           : " << Params.counter_max << "\n";
	cout << "Type of counters             : " << (Params.use_quake ? "Quake-compatibile\n" : "direct\n");
	if(Params.use_quake)
		cout << "Lowest quality value         : " << Params.lowest_quality << "\n";
	cout << "Both strands                 : " << (Params.both_strands ? "true\n" : "false\n");	
	cout << "RAM olny mode                : " << (Params.mem_mode ? "true\n" : "false\n");

	cout << "\n******* Stage 1 configuration: *******\n";
	cout << "\n";
	cout << "No. of bins                  : " << Params.n_bins << "\n";
	cout << "Bin part size                : " << Params.bin_part_size << "\n";
	cout << "Input buffer size            : " << Params.fastq_buffer_size << "\n";
	cout << "\n";

	cout << "No. of readers               : " << Params.n_readers << "\n";
	cout << "No. of splitters             : " << Params.n_splitters << "\n";
	cout << "\n";

	cout << "Max. mem. size               : " << setw(5) << (Params.max_mem_size / 1000000) << "MB\n";
	cout << "Max. mem. per storer         : " << setw(5) << (Params.max_mem_storer / 1000000) << "MB\n";
	cout << "Max. mem. for single package : " << setw(5) << (Params.max_mem_storer_pkg / 1000000) << "MB\n";
	cout << "\n";

	cout << "Max. mem. for PMM (bin parts): " << setw(5) << (Params.mem_tot_pmm_bins / 1000000) << "MB\n";
	cout << "Max. mem. for PMM (FASTQ)    : " << setw(5) << (Params.mem_tot_pmm_fastq / 1000000) << "MB\n";
	cout << "Max. mem. for PMM (reads)    : " << setw(5) << (Params.mem_tot_pmm_reads / 1000000) << "MB\n";

	cout << "\n";
}

//----------------------------------------------------------------------------------
// Show the settings of the KMC (in verbose mode only)
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::ShowSettingsStage2()
{
	if (!Params.verbose)
		return;

	cout << "\n******* Stage 2 configuration: *******\n";

	cout << "No. of sorters               : " << Params.n_sorters << "\n";
	cout << "No. of sort. threads         : ";
	for (uint32 i = 0; i < Params.n_omp_threads.size() - 1; ++i)
		cout << Params.n_omp_threads[i] << ", ";
	cout << Params.n_omp_threads.back() << "\n";

	cout << "\n";

	cout << "Max. mem. for 2nd stage      : " << setw(5) << (Params.max_mem_stage2 / 1000000) << "MB\n";
	cout << "\n";	
}
//----------------------------------------------------------------------------------
// Run the counter
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> bool CKMC<KMER_T, SIZE, QUAKE_MODE>::Process()
{
	int32 bin_id;
	CMemDiskFile *file;
	string name;
	uint64 size;
	uint64 n_rec;
	uint64 n_plus_x_recs;
	uint64 n_super_kmers;

	if (!initialized)
		return false;

	if (!AdjustMemoryLimits())
		return false;
	

	w1.startTimer();

	// Create monitors
	Queues.mm = new CMemoryMonitor(Params.max_mem_stage2);


	// Create queues
	Queues.input_files_queue = new CInputFilesQueue(Params.input_file_names);
	Queues.part_queue = new CPartQueue(Params.n_readers);
	Queues.bpq = new CBinPartQueue(Params.n_splitters);
	Queues.bd = new CBinDesc;
	Queues.bq = new CBinQueue(1);

	Queues.stats_part_queue = new CStatsPartQueue(Params.n_readers, STATS_FASTQ_SIZE);

	// Create memory manager
	Queues.pmm_bins = new CMemoryPool(Params.mem_tot_pmm_bins, Params.mem_part_pmm_bins);
	Queues.pmm_fastq = new CMemoryPool(Params.mem_tot_pmm_fastq, Params.mem_part_pmm_fastq);
	Queues.pmm_reads = new CMemoryPool(Params.mem_tot_pmm_reads, Params.mem_part_pmm_reads);
	Queues.pmm_stats = new CMemoryPool(Params.mem_tot_pmm_stats, Params.mem_part_pmm_stats);

	

	Queues.s_mapper = new CSignatureMapper(Queues.pmm_stats, Params.signature_len);
	
	// ***** Stage 0 *****
	w0.startTimer();
	w_stats_splitters.resize(Params.n_splitters);
	

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_stats_splitters[i] = new CWStatsSplitter<false>(Params, Queues);
		gr0_2.push_back(thread(std::ref(*w_stats_splitters[i])));
	}

	w_stats_fastqs.resize(Params.n_readers);
	
	for (int i = 0; i < Params.n_readers; ++i)
	{
		w_stats_fastqs[i] = new CWStatsFastqReader(Params, Queues);
		gr0_1.push_back(thread(std::ref(*w_stats_fastqs[i])));
	}
	for (auto p = gr0_1.begin(); p != gr0_1.end(); ++p)
		p->join();
	for (auto p = gr0_2.begin(); p != gr0_2.end(); ++p)
		p->join();


	uint32 *stats;
	Queues.pmm_stats->reserve(stats);
	fill_n(stats, (1 << Params.signature_len * 2) + 1, 0);


	for (int i = 0; i < Params.n_readers; ++i)
		delete w_stats_fastqs[i];

	for (int i = 0; i < Params.n_splitters; ++i)
	{
		w_stats_splitters[i]->GetStats(stats);			
		delete w_stats_splitters[i];
	}		

	delete Queues.stats_part_queue;
	Queues.stats_part_queue = NULL;
	delete Queues.input_files_queue;
	Queues.input_files_queue = new CInputFilesQueue(Params.input_file_names);

	heuristic_time.startTimer();
	Queues.s_mapper->Init(stats);
	heuristic_time.stopTimer();

	cout << "\n";
	
	w0.stopTimer();


	Queues.pmm_stats->free(stats);
	Queues.pmm_stats->release();
	delete Queues.pmm_stats;
	Queues.pmm_stats = NULL;

	Params.n_bins = Queues.s_mapper->get_max_bin_no() + 1;

	// ***** Stage 1 *****
	ShowSettingsStage1();

	w_splitters.resize(Params.n_splitters);

	for(int i = 0; i < Params.n_splitters; ++i)
	{
		w_splitters[i] = new CWSplitter<QUAKE_MODE>(Params, Queues);
		gr1_2.push_back(thread(std::ref(*w_splitters[i])));
	}

	w_storer = new CWKmerBinStorer(Params, Queues);
	gr1_3.push_back(thread(std::ref(*w_storer)));

	w_fastqs.resize(Params.n_readers);
	for(int i = 0; i < Params.n_readers; ++i)
	{
		w_fastqs[i] = new CWFastqReader(Params, Queues);
		gr1_1.push_back(thread(std::ref(*w_fastqs[i])));
	}

	for(auto p = gr1_1.begin(); p != gr1_1.end(); ++p)
		p->join();
	for(auto p = gr1_2.begin(); p != gr1_2.end(); ++p)
		p->join();

	Queues.pmm_fastq->release();
	Queues.pmm_reads->release();
	
	delete Queues.pmm_fastq;
	delete Queues.pmm_reads;

	for(auto p = gr1_3.begin(); p != gr1_3.end(); ++p)
		p->join();

	n_reads = 0;

	thread *release_thr_st1_1 = new thread([&]{
		for(int i = 0; i < Params.n_readers; ++i)
			delete w_fastqs[i];

		for(int i = 0; i < Params.n_splitters; ++i)
		{
			uint64 _n_reads;
			w_splitters[i]->GetTotal(_n_reads);
			n_reads += _n_reads;
			delete w_splitters[i];
		}

		delete w_storer;
	});

	thread *release_thr_st1_2 = new thread([&]{
		Queues.pmm_bins->release();
		delete Queues.pmm_bins;
	});


	release_thr_st1_1->join();
	release_thr_st1_2->join();

	delete release_thr_st1_1;
	delete release_thr_st1_2;


	w1.stopTimer();
	w2.startTimer();
	

	// ***** End of Stage 1 *****

	// Adjust RAM for 2nd stage
	// Calculate LUT size
	uint32 best_lut_prefix_len = 0;
	uint64 best_mem_amount = 1ull << 62;

	for (Params.lut_prefix_len = 2; Params.lut_prefix_len < 16; ++Params.lut_prefix_len)
	{
		uint32 suffix_len = Params.kmer_len - Params.lut_prefix_len;
		if (suffix_len % 4)
			continue;

		uint64 est_suf_mem = n_reads * suffix_len;
		uint64 lut_mem = Params.n_bins * (1ull << (2 * Params.lut_prefix_len)) * sizeof(uint64);

		if (est_suf_mem + lut_mem < best_mem_amount)
		{
			best_lut_prefix_len = Params.lut_prefix_len;
			best_mem_amount = est_suf_mem + lut_mem;
		}
	}

	Params.lut_prefix_len = best_lut_prefix_len;

#ifdef DEVELOP_MODE
	save_bins_stats(Queues, Params, sizeof(KMER_T), KMER_T::QUALITY_SIZE, n_reads);
#endif




	

	Queues.bd->reset_reading();
	vector<int64> bin_sizes;

	while((bin_id = Queues.bd->get_next_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);
		if (Params.max_x)
			bin_sizes.push_back(n_plus_x_recs * 2 * sizeof(KMER_T));			// estimation of RAM for sorting bins
		else
			bin_sizes.push_back(n_rec * 2 * sizeof(KMER_T));
	}
	
	sort(bin_sizes.begin(), bin_sizes.end(), greater<int64>());
	
	
	
	SetThreads2Stage(bin_sizes);
	AdjustMemoryLimitsStage2();

	Queues.kq = new CKmerQueue(Params.n_bins, Params.n_sorters);
	
	int64 stage2_size = 0;
	for (int i = 0; i < 4 * Params.n_sorters; ++i)
		stage2_size += bin_sizes[i];
	stage2_size = MAX(stage2_size, 16 << 20);
	Params.max_mem_stage2 = MIN(Params.max_mem_stage2, stage2_size);

	ShowSettingsStage2();
	
	// ***** Stage 2 *****
	Queues.bd->reset_reading();
	Queues.pmm_radix_buf  = new CMemoryPool(Params.mem_tot_pmm_radix_buf, Params.mem_part_pmm_radix_buf );
	if (!Params.use_quake && Params.both_strands)
		Queues.pmm_expand = new CMemoryPool(Params.mem_tot_pmm_epxand, Params.mem_part_pmm_epxand);
	else
		Queues.pmm_expand = NULL;
	Queues.memory_bins    = new CMemoryBins(Params.max_mem_stage2, Params.n_bins);
	if (Params.use_quake)
		Queues.pmm_prob = new CMemoryPool(Params.mem_tot_pmm_prob, Params.mem_part_pmm_prob);
	else
		Queues.pmm_prob = NULL;
	w_reader = new CWKmerBinReader<KMER_T, SIZE>(Params, Queues);
	gr2_1.push_back(thread(std::ref(*w_reader)));

	w_sorters.resize(Params.n_sorters);
	

	for(int i = 0; i < Params.n_sorters; ++i)
	{
		w_sorters[i] = new CWKmerBinSorter<KMER_T, SIZE>(Params, Queues, i);
		gr2_2.push_back(thread(std::ref(*w_sorters[i])));
	}

	w_completer = new CWKmerBinCompleter(Params, Queues);
	gr2_3.push_back(thread(std::ref(*w_completer)));

	for(auto p = gr2_1.begin(); p != gr2_1.end(); ++p)
		p->join();
	for(auto p = gr2_2.begin(); p != gr2_2.end(); ++p)
		p->join();


	for(auto p = gr2_3.begin(); p != gr2_3.end(); ++p)
		p->join();


	// ***** End of Stage 2 *****
	w_completer->GetTotal(n_unique, n_cutoff_min, n_cutoff_max, n_total);

	thread *release_thr_st2_1 = new thread([&]{
		delete Queues.mm;
		if (Queues.pmm_expand)
		{
			Queues.pmm_expand->release();
			delete Queues.pmm_expand;
		}
		Queues.pmm_radix_buf->release();
		Queues.memory_bins->release();
		delete Queues.pmm_radix_buf;
		delete Queues.memory_bins;
	});

	
	uint64 stat_n_plus_x_recs, stat_n_recs, stat_n_recs_tmp, stat_n_plus_x_recs_tmp;
	stat_n_plus_x_recs = stat_n_recs = stat_n_recs_tmp = stat_n_plus_x_recs_tmp = 0;
	thread *release_thr_st2_2 = new thread([&]{
		
		delete w_reader;
		for(int i = 0; i < Params.n_sorters; ++i)
		{
			w_sorters[i]->GetDebugStats(stat_n_recs_tmp, stat_n_plus_x_recs_tmp);
			stat_n_plus_x_recs += stat_n_plus_x_recs_tmp;
			stat_n_recs += stat_n_recs_tmp;
			delete w_sorters[i];
		}
		delete w_completer;

		delete Queues.input_files_queue;
		delete Queues.bq;
		delete Queues.part_queue;
		delete Queues.bpq;
		delete Queues.kq;
	});

	


	// ***** Removing temporary files *****

	tmp_size = 0;
	n_total_super_kmers = 0;
	Queues.bd->reset_reading();
	while((bin_id = Queues.bd->get_next_bin()) >= 0)
	{
		Queues.bd->read(bin_id, file, name, size, n_rec, n_plus_x_recs, n_super_kmers);		
#ifndef DEVELOP_MODE
		boost::filesystem::remove(boost::filesystem::path(name));
#endif // DEVELOP_MODE
		tmp_size += size;
		n_total_super_kmers += n_super_kmers;
	}
	delete Queues.bd;

	release_thr_st2_1->join();
	release_thr_st2_2->join();

	//KMC_2
	//cout << "n_recs: " << stat_n_recs << ", n_plus_x_recs: " << stat_n_plus_x_recs << "(" << (float)(stat_n_plus_x_recs) / stat_n_recs << ")\n";

	delete release_thr_st2_1;
	delete release_thr_st2_2;
	delete Queues.s_mapper;
	w2.stopTimer();

	return true;
}

//----------------------------------------------------------------------------------
// Return statistics
template <typename KMER_T, unsigned SIZE, bool QUAKE_MODE> void CKMC<KMER_T, SIZE, QUAKE_MODE>::GetStats(double &time1,
	double &time2, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64& _n_total_super_kmers)
{
	time1 = w1.getElapsedTime();
	time2 = w2.getElapsedTime();

	_n_unique	  = n_unique;
	_n_cutoff_min = n_cutoff_min;
	_n_cutoff_max = n_cutoff_max;
	_n_total      = n_total;
	_n_reads      = n_reads;
	_tmp_size     = tmp_size;
	_n_total_super_kmers = n_total_super_kmers;
}

#endif

// ***** EOF
