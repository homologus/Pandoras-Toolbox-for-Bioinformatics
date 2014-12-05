#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>

#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)

int max_memory; // the most memory one should alloc at any time, in MB

int order = 0; // deblooming order; 0 = debloom everything; 1 = don't debloom 1-node tips (experimental, untested, shouldn't work);// (made extern int in Traversal.h)
  
#include "Bank.h"
#include "Hash16.h"
#include "Set.h"
#include "Pool.h"
#include "Bloom.h"
#include "Debloom.h"
#include "Utils.h"
#include "SortingCount.h"
#include "Terminator.h"
#include "Kmer.h"
#include "Traversal.h"
#include "rvalues.h" // for 4bloom


int64_t genome_size;
Bloom * bloo1;


inline void assemble()
{

    //////-------------------------------------------------------------------------------------------
    fprintf (stderr,"______________________________________________________ \n");
    fprintf (stderr,"___________ Assemble from bloom filter _______________ \n");
    fprintf (stderr,"______________________________________________________ \n\n");

    //////-------------------------------------------------------------------------------------------


    long long len_left = 0;
    long long len_right = 0;
    long long contig_len =0;
    long long maxlen=10000000;

    char *left_traversal  = (char *) malloc(maxlen*sizeof(char));
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    char *contig          = (char *) malloc(2*(maxlen+sizeKmer)*sizeof(char));
    kmer_type kmer;

    long long nbContig =0;
    long long nbSmallContig =0;
    long long totalnt=0;
    long long max_contig_len=0;
    long long mlenleft=0,mlenright=0;
    int64_t NbBranchingKmer=0;
    char kmer_seq[sizeKmer+1];
    FILE * file_assembly = fopen(return_file_name(assembly_file),"w+");

    BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);

    STARTWALL(assembly);

    char *assemble_only_one_region = NULL; // debugging, set to a ASCII kmer to activate, NULL to desactivate
    bool LOAD_BRANCHING_KMERS=false; // debugging
    bool DUMP_BRANCHING_KMERS=false;
   
    BranchingTerminator *terminator;

    if (LOAD_BRANCHING_KMERS)
    {
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),false);
        terminator = new BranchingTerminator(BranchingKmers,SolidKmers, bloo1,false_positives);
        BranchingKmers->close();
    }
    else
        terminator = new BranchingTerminator(SolidKmers,genome_size, bloo1,false_positives);

    if (DUMP_BRANCHING_KMERS)
    {
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),true);
        terminator->dump_branching_kmers(BranchingKmers);
        BranchingKmers->close();
    }

#ifdef UNITIG
    SimplePathsTraversal *traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    fprintf (stderr,"_________________Assembling in Unitig mode ..._____________________ \n\n");
#else
    MonumentTraversal *traversal = new MonumentTraversal(bloo1,false_positives,terminator);
#endif
    //RandomBranchingTraversal *traversal = new RandomBranchingTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(maxlen);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    
    while (terminator->next(&kmer))
    {
        // keep looping while a starting kmer is available from this kmer
		// everything will be marked during the traversal()'s
		kmer_type starting_kmer;
#ifdef UNITIG
        while (traversal->get_new_starting_node_improved(kmer,starting_kmer))
#else
        while (traversal->find_starting_kmer(kmer,starting_kmer))
#endif
		{
		    code2seq(starting_kmer,kmer_seq); // convert starting kmer to nucleotide seq
            traversal->revert_stats(); // set stats from the last commit (discard stats from find_starting_kmer / small contigs)

            if (assemble_only_one_region != NULL)
            {
                kmer_type dummy;
                starting_kmer = extractKmerFromRead(assemble_only_one_region,0,&kmer,&dummy,false);
            }

            // right extension
            len_right = traversal->traverse(starting_kmer,right_traversal,0);
            mlenright= max(len_right,mlenright);

            // left extension, is equivalent to right extension of the revcomp
            len_left = traversal->traverse(starting_kmer,left_traversal,1);
            mlenleft= max(len_left,mlenleft);

            // form the contig
            revcomp_sequence(left_traversal,len_left);
            strcpy(contig,left_traversal); // contig = revcomp(left_traversal)
	        strcat(contig,kmer_seq);//               + starting_kmer
            strcat(contig,right_traversal);//           + right_traversal

            contig_len=len_left+len_right+sizeKmer;

            // save the contig
            if(contig_len >= MIN_CONTIG_SIZE)
            {
                max_contig_len = max(max_contig_len,contig_len);
                fprintf(file_assembly,">%lli__len__%lli \n",nbContig,contig_len);
                fprintf(file_assembly,"%s\n",contig);
                nbContig++;
                totalnt+=contig_len;
                traversal->commit_stats();
            }
            else
            {
                traversal->revert_stats();
                nbSmallContig++;
            }
            if (assemble_only_one_region != NULL)
                break;
        }
    
        NbBranchingKmer++;
        if ((NbBranchingKmer%300)==0) fprintf (stderr,"%cLooping through branching kmer n° %lld / %lld  total nt   %lld   ",13,(long long int) NbBranchingKmer,(long long int) terminator->nb_branching_kmers, (long long int)totalnt );

        if (nbContig > 0 && assemble_only_one_region != NULL)
            break;

    }
    fclose(file_assembly);

    fprintf (stderr,"\n Total nt assembled  %lli  nbContig %lli\n",totalnt,nbContig);
    fprintf (stderr," Max contig len  %lli (debug: max len left %lli, max len right %lli)\n",max_contig_len,mlenleft,mlenright);
    fprintf (stderr,"\n Debug traversal stats: %ld ends of contigs (%lld unsaved small contigs), among them:\n",traversal->final_stats.ended_traversals,nbSmallContig);
    fprintf (stderr," %ld couldn't validate consensuses\n",traversal->final_stats.couldnt_validate_consensuses);
    fprintf (stderr," %ld large bubble breadth, %ld large bubble depth, %ld marked kmer, %ld no extension\n",traversal->final_stats.couldnt_traverse_bubble_breadth,traversal->final_stats.couldnt_traverse_bubble_depth,traversal->final_stats.couldnt_because_marked_kmer,traversal->final_stats.couldnt_find_extension);
    fprintf (stderr," %ld in-branchin large depth, %ld in-branching large breadth, %ld in-branching other\n",traversal->final_stats.couldnt_inbranching_depth,traversal->final_stats.couldnt_inbranching_breadth,traversal->final_stats.couldnt_inbranching_other);
    
    STOPWALL(assembly,"Assembly");

    free(left_traversal);
    free(right_traversal);
    free(contig);
    SolidKmers->close();
}

int main(int argc, char *argv[])
{
    
    if(argc <  6)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size min_abundance estimated_genome_size prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_genome_size is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    bool FOUR_BLOOM_VERSION = true;

     // shortcuts to go directly to assembly using serialized bloom and serialized hash
    int START_FROM_SOLID_KMERS=0; // if = 0, construct the fasta file of solid kmers, if = 1, start directly from that file 
    int LOAD_FALSE_POSITIVE_KMERS=0; // if = 0, construct the fasta file of false positive kmers (debloom), if = 1, load that file into the hashtable
    int NO_FALSE_POSITIVES_AT_ALL=0; // if = 0, normal behavior, if = 1, don't load false positives (will be a probabilistic de bruijn graph)
    int max_disk_space = 0;// let dsk decide
    for (int n_a = 6; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"--original") == 0)
    	    FOUR_BLOOM_VERSION = false;

        if (strcmp(argv[n_a],"--dont-count")==0)
            START_FROM_SOLID_KMERS = 1;

        if (strcmp(argv[n_a],"--dont-debloom")==0)
            LOAD_FALSE_POSITIVE_KMERS = 1;

        if (strcmp(argv[n_a],"--just-assemble")==0)
        {
            START_FROM_SOLID_KMERS = 1;
            LOAD_FALSE_POSITIVE_KMERS = 1;
        }

        if (strcmp(argv[n_a],"--titus-mode")==0)
            NO_FALSE_POSITIVES_AT_ALL = 1;
        
        
        if (strcmp(argv[n_a],"-d")==0)
            max_disk_space = atoi(argv[n_a+1]);
        
        
        if (strcmp(argv[n_a],"-maxc")==0)
	    max_couv = atoi(argv[n_a+1]);
        
        if (strcmp(argv[n_a],"--le-changement")==0)
            {printf("c'est maintenant!\n");exit(0);}
    }


    // kmer size
    sizeKmer=27; // let's make it even for now, because i havnt thought of how to handle palindromes (dont want to stop on them)
    if(argc >=  3)
    {
        sizeKmer = atoi(argv[2]);
        if (sizeKmer%2==0)
        {
            sizeKmer-=1;
            printf("Need odd kmer size to avoid palindromes. I've set kmer size to %d.\n",sizeKmer);
        }
        if (sizeKmer>((int)sizeof(kmer_type)*4))
        {
            printf("Max kmer size on this compiled version is %lu\n",sizeof(kmer_type)*4);
            exit(1);
        }
    }

    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    double lg2 = log(2);
   
    if (sizeKmer > 128)
    {
        FOUR_BLOOM_VERSION = false;
        printf("Reverted to single Bloom filter implementation for k>128\n");
    }

    if (!FOUR_BLOOM_VERSION) 
      NBITS_PER_KMER = log(16*sizeKmer*(lg2*lg2))/(lg2*lg2); // needed to process argv[5]
    else 
      NBITS_PER_KMER = rvalues[sizeKmer][1];

    // solidity 
    nks =NNKS;
    if(argc >=  4)
    {
        nks = atoi(argv[3]);
        if (nks==0) nks=1; // min abundance can't be 0
    }


   if(argc >=  5)
    {
       genome_size  = atoll(argv[4]);
      // int estimated_bloom_size = max( (int)ceilf(log2f(genome_size * NBITS_PER_KMER )), 1);
        uint64_t estimated_bloom_size = (uint64_t) (genome_size * NBITS_PER_KMER);

       uint64_t estimated_nb_FP =  (uint64_t)(genome_size * 4 * powf(0.6,11)); // just indicative
    
       //max_memory = max( (1LL << estimated_bloom_size)/8LL /1024LL/1024LL, 1LL );
        max_memory =  max((int64_t) estimated_bloom_size/8LL /1024LL/1024LL,1LL);

      printf("estimated values: nbits Bloom %lli, nb FP %lld, max memory %i MB\n",estimated_bloom_size,estimated_nb_FP,max_memory);

    }

    // output prefix
    if(argc >=  6)
    {
        strcpy(prefix,argv[5]);
    }

   


    fprintf (stderr,"taille cell %lu \n", sizeof(cell<kmer_type>));

    STARTWALL(0);

    Bank *Reads = new Bank(argv[1]);
    
    // counter kmers, write solid kmers to disk
    if (!START_FROM_SOLID_KMERS)
    {
        int verbose = 0;
        bool write_count = false;
        bool skip_binary_conversion = false;

        sorting_count(Reads,prefix,max_memory,max_disk_space,write_count,verbose, skip_binary_conversion);
    }

    // debloom, write false positives to disk, insert them into false_positives
    if (! LOAD_FALSE_POSITIVE_KMERS)
    {
        debloom(order, max_memory);
    }
    
    bloo1 = bloom_create_bloo1((BloomCpt *)NULL, false);

    if (! NO_FALSE_POSITIVES_AT_ALL)
    {
        // load false positives from disk into false_positives
        if (!FOUR_BLOOM_VERSION) 
            false_positives = load_false_positives();
	else
	    false_positives = load_false_positives_cascading4();
    }
    else
    {
        // titus mode: no FP's
        false_positives = dummy_false_positives();
    }

    //  return 1;
    assemble(); 

    STOPWALL(0,"Total");

    delete Reads;
    return 0;
}


