//
//  Terminator.cpp
//

#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif 

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm> // for max

#include "Terminator.h"
#include "Traversal.h" // for extensions() 

using namespace::std;
bool Terminator::verbose = true;

// common terminator functions (actually, more like Kmer64 + bloom + debloom)

// we define a structure 2x4 bits which indicates which nucleotides correspond to in- and out-branching of a kmer
// let's see an example (bidirected debruijn graph):
// 
//     ACT [G]  -\/->   [C] CAG [T]  <--  [C] CCA
//     ---       X          ---               ---
//[C]  AGT     -/\- >   [A] CTG [G]   -->     TGG [G]
//
// structure for the node CAG/CTG: (forward:)0010(reverse:)0001 (the order is ACTG); T forward (to the right) and G reverse (to the right)
unsigned char Terminator::branching_structure(kmer_type graine)
{
    assert(graine<=revcomp(graine));
    kmer_type new_graine;
    unsigned char result = 0;
    int nt;
    int strand;

    for(nt=0; nt<4; nt++) 
    {
        // forward right extensions 
        strand=0;
        new_graine = next_kmer(graine,nt,&strand);
        if(bloom_solid_kmers->contains(new_graine) && !debloom->contains(new_graine)){
            result|=1<<nt;
        }

        // reverse right extensions
        strand=1;
        new_graine = next_kmer(graine,nt,&strand);
        if(bloom_solid_kmers->contains(new_graine) && !debloom->contains(new_graine)){
            result|=1<<(nt+4);
        }
    }
    return result;
}

// determines if a kmer is branching or not
bool Terminator::is_branching(kmer_type graine)
{
    assert(graine<=revcomp(graine));

    // method specific to order=0 (remember that order>0 is never used)
    if (order == 0)
    {
        // cannot really be optimized, because most kmers will be non-branching, hence computing branching_structure() takes optimal time
        int nb_forward_links = 0, nb_reverse_links = 0;
        int i;
        unsigned char branching = branching_structure(graine);

        for (i=0;i<4;i++)
            nb_forward_links += (branching>>i)&1;

        for (i=4;i<8;i++)
            nb_reverse_links += (branching>>i)&1;

        return !(nb_forward_links == 1 && nb_reverse_links == 1);
    }
    else
    {
        if (order==1 && is_tip(graine,bloom_solid_kmers,debloom))// algo fixme-order>0: order=1, tips are ignored as branching kmers. should find a more general code for order>1
            return false;

        // any order>=0: check that this node really yields a true branching, as opposed to branching yielding tips
        for (int strand=0;strand<2;strand++)
        {
            int osef;
            int nb_extensions = traversal_extensions(graine,strand,osef, bloom_solid_kmers, debloom);
            if (nb_extensions != 1)
                return true;
        }
    }
    return false;
}

// [branching kmers]-based terminator
// most kmers should have 1 in-neighbor and 1 out-neighbor
// this terminator indexes all the other kmers (ie., branching or dead-end kmers)
// it will miss circular regions tho..

// mark a kmer (the kmer has to be a branching kmer or a neighbor of one, else no effect)
void BranchingTerminator::mark(kmer_type graine)
{
    assert(graine<=revcomp(graine));
    bool could_mark = false;

    // if it is a branching kmer, mark it directly (it may have no branching neighbor)
    if (is_indexed(graine))
    {
        set_value_t val;
        branching_kmers->get(graine,&val);
        branching_kmers->set(graine,val|(1<<8));
        could_mark = true;
    }

    // enumerate all the neighbors
    for (int strand=0; strand < 2; strand++)
    {
        for(int nt=0; nt<4; nt++) 
        {
            // test if the neighbor is a branching kmer
            int neighbor_strand = strand;
            kmer_type neighbor = next_kmer(graine,nt,&neighbor_strand);
            if( (!bloom_solid_kmers->contains(neighbor)) || debloom->contains(neighbor) )
                continue;
            if (!is_indexed(neighbor))
                continue;

            // mark the kmer in that neighbor
            int revStrand ;
            int revNT;
            revStrand = 1-neighbor_strand;
            if (strand == 0)
                revNT = revcomp_int(code2nucleotide(graine,0));
            else
                revNT = code2nucleotide(graine,sizeKmer-1);

            mark(neighbor,revNT,revStrand);

            could_mark = true;
        }
    }
    if (could_mark)
        assert(is_marked(graine));
}

// record info into a branching kmer
void BranchingTerminator::mark(kmer_type graine, char nt, int strand)
{
    assert(nt<4);
    assert(strand<2);
    assert(graine<=revcomp(graine));

    // BranchingTerminator ignores non-branching kmers
    if (!is_indexed(graine))
        return;

    //int val = 0;
    set_value_t val=0;
    branching_kmers->get(graine,&val);

    // set a 1 at the right NT & strand position
    if (strand==0)
        val|=1<<(nt);
    else
        val|=1<<(nt+4);

    //    printf ("mark,  graine = %llx val: %x branching structure: %x\n",graine,val,branching_structure(graine));
    
    branching_kmers->set(graine,val); //was insert for Hash16

    assert(is_marked(graine,nt,strand));
}

// test if a kmer is marked, providing that the kmer is a branching kmer or a neighbor of one (else, considered unmarked)
bool BranchingTerminator::is_marked(kmer_type graine)
{
    assert(graine<=revcomp(graine));

    // if it is a branching kmer, read marking directly (it may have no branching neighbor)
    if (is_indexed(graine))
        return is_marked_branching(graine);

    // enumerate all the neighbors
    for (int strand=0; strand < 2; strand++)
    {
        for(int nt=0; nt<4; nt++) 
        {
            // test if the neighbor is a branching kmer
            int neighbor_strand = strand;
            kmer_type neighbor = next_kmer(graine,nt,&neighbor_strand);
            if( (!bloom_solid_kmers->contains(neighbor)) || debloom->contains(neighbor) )
                continue;
            if ( !is_indexed(neighbor) )
                continue;

            // test the kmer w.r.t that neighbor
            int revStrand ;
            int revNT;
            revStrand = 1-neighbor_strand;
            if (strand == 0)
                revNT = revcomp_int(code2nucleotide(graine,0));
            else
                revNT = code2nucleotide(graine,sizeKmer-1);

            if ( is_marked(neighbor,revNT,revStrand) ) 
                return true;
        }
    }
    return false;
}

// test if a branching kmer is marked (using the special marker for branching kmers only)
bool BranchingTerminator::is_marked_branching(kmer_type graine)
{
    assert(graine<=revcomp(graine));
    assert(is_branching(graine));
    assert(is_indexed(graine));
    
    set_value_t val;
    branching_kmers->get(graine,&val);
    return (val&(1<<8)) != 0;
}

// that function returns false for non-indexed kmers
bool BranchingTerminator::is_marked(kmer_type graine, char nt, int strand)
{
    assert(nt<4);
    assert(strand<2);
    assert(graine<=revcomp(graine));

    set_value_t val = 0;
    int is_present = branching_kmers->get(graine,&val);

    if (!is_present)
        return false;

    //printf ("is_marked,  graine = %llx val: %x branching structure: %x\n",graine,val,branching_structure(graine));

    int extension_nucleotide_marked;
    if (strand==0)
        extension_nucleotide_marked = (val>>nt)&1;
    else
        extension_nucleotide_marked = (val>>(nt+4))&1;

    return  extension_nucleotide_marked == 1;
}

bool BranchingTerminator::is_indexed(kmer_type graine)
{
  return branching_kmers->contains(graine);
}

BranchingTerminator::BranchingTerminator(BinaryBank *given_SolidKmers, uint64_t genome_size, Bloom *given_bloom, Set *given_debloom) : Terminator(given_SolidKmers,given_bloom,given_debloom)
{
    // estimate, from the first million of kmers, the number of branching kmers, extrapolating given the estimated genome size
    // TODO: erwan noticed that this code isn't useful anymore with AssocSet, feel free to remove it sometimes
    uint64_t nb_extrapolation = 3000000;
    SolidKmers->rewind_all();
    uint64_t nb_kmers = 0;
    nb_branching_kmers = 0;
    kmer_type kmer;
    uint64_t previous_estimated_nb_branching_kmers, estimated_nb_branching_kmers;
    while (SolidKmers->read_element(&kmer))
    {
        if (is_branching(kmer))
            nb_branching_kmers++;

        if ((nb_branching_kmers%1000)==0 && nb_branching_kmers>0)
        {
            previous_estimated_nb_branching_kmers = estimated_nb_branching_kmers;
            estimated_nb_branching_kmers = (uint64_t)((1.0*nb_branching_kmers)/nb_kmers * genome_size);
            // minor todo: stop when previous_.. - estimated < threshold (pourquoi pas = 10% estimated)
            fprintf (stderr,"%cExtrapolating the number of branching kmers from the first %dM kmers: %lld",13,(int)ceilf(nb_extrapolation/1024.0/1024.0),estimated_nb_branching_kmers);
        }

        if (nb_kmers++ == nb_extrapolation)
            break;
    }
    estimated_nb_branching_kmers = (uint64_t)((1.0*nb_branching_kmers)/nb_kmers * genome_size); // final estimation
    int estimated_NBITS_TERMINATOR = max( (int)ceilf(log2f(estimated_nb_branching_kmers)), 1);
    fprintf (stderr,"\n");

    // call Hash16 constructor
    // branching_kmers = new Hash16(estimated_NBITS_TERMINATOR);
    branching_kmers = new AssocSet();

    // index, once and for all, all the branching solid kmers
    SolidKmers->rewind_all();
    nb_branching_kmers = 0;
    uint64_t nb_solid_kmers = 0;
    while (SolidKmers->read_element(&kmer))
    {
        if (is_branching(kmer))
        {
	  // branching_kmers->insert(kmer,0);
            branching_kmers->insert(kmer);
            nb_branching_kmers++;
        }

        nb_solid_kmers++;
        if ((nb_branching_kmers%500)==0) fprintf (stderr,"%cIndexing branching kmers %lld / ~%lld",13,nb_branching_kmers,estimated_nb_branching_kmers);
    }

    if (nb_branching_kmers == 0)
        printf("\n**** Warning\n\nNo branching kmers were found in this dataset (it is either empty or a tiny circular genome) - Minia will not assemble anything.\n\n****\n\n");

	branching_kmers->finalize();
//	branching_kmers->print_total_size();
//    fprintf (stderr,"\n\nAllocated memory for marking: %lld branching kmers x (%lu+%lu) B \n",nb_branching_kmers,sizeof(kmer_type),sizeof(set_value_t));
//    fprintf (stderr," actual implementation:  (sizeof(kmer_type) = %lu B) + (sizeof(set_value_t) = %lu B) per entry:  %.2f bits / solid kmer\n",sizeof(kmer_type),sizeof(set_value_t),(nb_branching_kmers*((sizeof(kmer_type)+sizeof(set_value_t))*8.0))/nb_solid_kmers);

    // init branching_kmers iterator for what happens next
     branching_kmers->start_iterator(); 
}

// constructor that simply loads a dump of branching kmers
BranchingTerminator::BranchingTerminator(BinaryBank *branchingKmers, BinaryBank *given_SolidKmers, Bloom *given_bloom, Set *given_debloom) : Terminator(given_SolidKmers,given_bloom,given_debloom)
{
    nb_branching_kmers = branchingKmers->nb_elements();
    int NBITS_TERMINATOR = max( (int)ceilf(log2f(nb_branching_kmers)), 1);

    // call Hash16 constructor
    // branching_kmers = new Hash16(NBITS_TERMINATOR);
    branching_kmers = new AssocSet();

    // load branching kmers
    branchingKmers->rewind_all();
    kmer_type kmer;
    while (branchingKmers->read_element(&kmer))
      branching_kmers->insert(kmer); //,0 for Hash16

	branching_kmers->finalize();

    if (verbose)
        fprintf (stderr,"\nLoaded %lld branching kmers x %lu B =  %.1f MB\n",nb_branching_kmers,sizeof(cell<kmer_type>),((1<<NBITS_TERMINATOR)*sizeof(cell<kmer_type>)*1.0)/1024.0/1024.0);

    // init branching_kmers iterator for what happens next
    branching_kmers->start_iterator(); 
}

BranchingTerminator::~BranchingTerminator()
{
    delete branching_kmers;
}

bool BranchingTerminator::next(kmer_type *kmer)
{   
    if (branching_kmers->next_iterator())
    {
      //*kmer =  branching_kmers->iterator.cell_ptr->graine; //for Hash16
      *kmer =  *(branching_kmers->iterator);

        return true;
    }
    return false;
}

void BranchingTerminator::dump_branching_kmers(BinaryBank *BranchingKmers)
{
    // init branching_kmers iterator for what happens next
    branching_kmers->start_iterator(); 
    
    while (branching_kmers->next_iterator())
    {
      //kmer_type kmer =  branching_kmers->iterator.cell_ptr->graine;//for Hash16
      kmer_type kmer = *( branching_kmers->iterator);

        BranchingKmers->write_element(&kmer);
    }
    fprintf (stderr,"Dumped branching kmers\n");
}


// reset all marking information: everything is now unmarked
void BranchingTerminator::reset()
{
    branching_kmers->clear(); 
}

//--------------------------------------------------------


// bloom-based terminator (untested, totally unused)
BloomTerminator::BloomTerminator(int tai_Bloom)
{
    bloo2 = new Bloom(tai_Bloom); // to mark kmers already used for assembly  ... taille a changer
}

void BloomTerminator::mark(kmer_type graine, char nt, int strand)
{
    bloo2->add(graine); 
}

bool BloomTerminator::is_marked(kmer_type graine, char nt, int strand)
{
    return bloo2->contains(graine);
}

bool BloomTerminator::is_fully_marked(kmer_type graine)
{
    return is_marked(graine,0,0);
}

bool BloomTerminator::next(kmer_type *kmer)
{    
    return SolidKmers->read_element(kmer);
}

