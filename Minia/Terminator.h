//
//  Terminator.h

#ifndef Terminator_h
#define Terminator_h
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <cmath> // for log2f

#include "Bloom.h"
#include "Kmer.h"
#include "Set.h"
#include "Bank.h"

class Terminator{

protected:

    BinaryBank *SolidKmers;
    Bloom *bloom_solid_kmers;  
    Set *debloom;

    unsigned char branching_structure(kmer_type graine);
public:
  static bool verbose;
  virtual void mark(kmer_type graine, char nt, int strand){ return; };
  virtual bool is_marked(kmer_type graine, char nt, int strand){ return (1); };
  virtual void mark(kmer_type graine){ return; };
  virtual bool is_marked(kmer_type graine){ return (1); };
  virtual bool is_marked_branching(kmer_type graine) {return 1;};
  virtual void reset() {return;};
  bool is_branching(kmer_type graine);
  bool next(kmer_type *kmer);

    Terminator(BinaryBank *given_SolidKmers, Bloom *given_bloom, Set *given_debloom) : SolidKmers(given_SolidKmers), bloom_solid_kmers(given_bloom), debloom(given_debloom) { }
};


class BloomTerminator {
protected:
    Bloom * bloo2;
    BinaryBank *SolidKmers;
public:  // is there a way to not repeat the declaration of Terminator functions?
  bool is_fully_marked(kmer_type graine);
  void mark(kmer_type graine, char nt, int strand);
  bool is_marked(kmer_type graine, char nt, int strand);
  bool next(kmer_type *kmer);
 
    BloomTerminator(int tai_Bloom);
    ~BloomTerminator();
};

class BranchingTerminator:  public Terminator{
  bool is_indexed(kmer_type graine);
  int genome_size;
  // Hash16 *branching_kmers;
  AssocSet *branching_kmers;
public:
   bool is_fully_marked(kmer_type graine);
  void mark(kmer_type graine, char nt, int strand);
  bool is_marked(kmer_type graine, char nt, int strand);
  void mark(kmer_type graine);
  bool is_marked(kmer_type graine);
  bool is_marked_branching(kmer_type graine);
  bool next(kmer_type *kmer);
    int64_t nb_branching_kmers;
    void dump_branching_kmers(BinaryBank *BranchingKmers);
    void reset();

    BranchingTerminator(BinaryBank *given_SolidKmers, uint64_t genome_size, Bloom *given_bloom, Set *given_debloom);
    BranchingTerminator(BinaryBank *branchingKmers, BinaryBank *given_SolidKmers, Bloom *given_bloom, Set *given_debloom);
    ~BranchingTerminator();
};




#endif
