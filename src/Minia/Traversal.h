#ifndef Traversal_h
#define Traversal_h

#include <set>
#include <queue>
#include <string>
#include <cmath> // only for pow()

#include "Kmer.h"
#include "Terminator.h"
#include "Bloom.h"
#include "Set.h"
#include "Utils.h" // for needleman_wunch

using namespace std;

// types using in advanced traversal functions
struct kmer_strand_nt { 
    kmer_type kmer;
    int strand;
    int nt;
    kmer_strand_nt(kmer_type kmer, int strand, int nt)  : kmer(kmer), strand(strand), nt(nt) {}
    bool operator<(const kmer_strand_nt &other) const {
        // need to define a strict weak ordering
        if (kmer != other.kmer)
            return (kmer < other.kmer);
        if  (strand != other.strand)
            return (strand < other.strand);
        return (nt < other.nt); }
};

// in traversals, a node is just a tuple (kmer, strand, nt initially chosen to reach that node)
// so it just happens that it's the same type as kmer_strand_nt
typedef kmer_strand_nt node; 
typedef queue<node> queue_nodes;

// some stats
struct TraversalStats
{
    long ended_traversals;
    long couldnt_find_all_consensuses;
    long couldnt_validate_consensuses;
    long couldnt_traverse_bubble_breadth;
    long couldnt_traverse_bubble_depth;
    long couldnt_because_marked_kmer;
    long couldnt_inbranching_depth;
    long couldnt_inbranching_breadth;
    long couldnt_inbranching_other;
    long couldnt_find_extension;
};

// semi-abstract class. implements traverse but not avance 
class Traversal{

protected:
    Bloom *bloom;
    Set *debloom;
    Terminator *terminator;
    int maxlen;
    int max_depth;
    int max_breadth;

    virtual char avance(kmer_type graine, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer) = 0;

    void mark_extensions(set<kmer_type> *extensions_to_mark);

public:
    Traversal(Bloom *given_bloom, Set *given_debloom, Terminator *given_terminator) : bloom(given_bloom), debloom(given_debloom), terminator(given_terminator), maxlen(1000000),max_depth(500),max_breadth(20), stats(TraversalStats()), final_stats(TraversalStats()) { }
    ~Traversal();
    
    void set_maxlen(int);
    void set_max_depth(int);
    void set_max_breadth(int);
    int traverse(kmer_type starting_kmer, char* resulting_sequence, int current_strand, kmer_type previous_kmer = 0);

    // n-order extension function, to ignore tips
    int extensions(kmer_type kmer, int strand, int &nt);
    
    // useful atomic avance functions
    int simple_paths_avance(kmer_type graine, int current_strand, bool first_extension, char * newNT);
    char random_unmarked_avance(kmer_type graine, int current_strand, bool first_extension, char * newNT);
	
	// high level starting kmer selection
    bool get_new_starting_node(kmer_type branching_kmer, kmer_type &starting_kmer);
    bool get_new_starting_node_improved(kmer_type branching_kmer, kmer_type &starting_kmer);
	bool find_starting_kmer_inside_simple_path(kmer_type kmer, kmer_type &starting_kmer); // now unused

    vector<pair<int, int> > bubbles_positions; // record the start/end positions of traversed bubbles (only from the latest traverse() call)

    TraversalStats final_stats, stats;
    void commit_stats(); // save current stats into final stats
    void revert_stats(); // discard changes in stats (because discarded contig)
};

// n-order extension function, to ignore tips
extern int order; // declared and initialized in assemb.cpp
int traversal_extensions(kmer_type kmer, int strand, int &nt, Bloom *bloom_solid_kmers, Set *debloom);
bool is_tip(kmer_type kmer, Bloom *bloom_solid_kmers, Set *debloom);

class RandomBranchingTraversal: public Traversal
{
protected:
    char avance(kmer_type graine, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer);

public:
    RandomBranchingTraversal(Bloom *given_bloom, Set *given_debloom, Terminator *given_terminator) : Traversal(given_bloom,given_debloom,given_terminator) {}
    
};


class SimplePathsTraversal: public Traversal
{
    char avance(kmer_type graine, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer);
public:
    SimplePathsTraversal(Bloom *given_bloom, Set *given_debloom, Terminator *given_terminator) : Traversal(given_bloom,given_debloom,given_terminator) {}

};

// auxiliary class that is used by MonumentTraversal and deblooming 
class Frontline
{
    kmer_type starting_kmer; int starting_strand;  Bloom *bloom; Set *debloom; Terminator *terminator; set<kmer_type> *all_involved_extensions; kmer_type previous_kmer;
    queue_nodes frontline;
    // set<node, CompareWithoutNT> already_frontlined;
    set<kmer_type> already_frontlined; // making it simpler now
public:
    bool check_in_branching;
    int depth;
    Frontline(kmer_type starting_kmer, int starting_strand, Bloom *bloom, Set *debloom, Terminator *terminator, set<kmer_type> *all_involved_extensions, kmer_type previous_kmer = 0, bool check_in_branching = true); 
    bool go_next_depth();
    int size();
    node front();
    bool check_inbranching(kmer_type from_kmer, int current_strand);
    enum reason 
    {
        NONE,
        ALREADY_FRONTLINED,
        IN_BRANCHING_DEPTH,
        IN_BRANCHING_BREADTH,
        IN_BRANCHING_OTHER,
        MARKED
    };
    reason stopped_reason;

};

class MonumentTraversal: public Traversal
{
    //int max_length_deadend = 150; // replaced by sizeKmer+1
    static const int consensuses_identity = 90; // traversing bubble if paths are all pair-wise identical by > 90%
    char avance(kmer_type graine, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer);

    set<string> all_consensuses_between(kmer_type start_kmer, int start_strand, kmer_type end_kmer, int end_strand, int traversal_depth, set<kmer_type> used_kmers, string current_consensus, bool &success);
    set<string> all_consensuses_between(kmer_type start_kmer, int start_strand, kmer_type end_kmer, int end_node, int traversal_depth, bool &success);
    int find_end_of_branching(kmer_type starting_kmer, int starting_strand, kmer_type &end_kmer, int &end_strand, kmer_type previous_kmer, set<kmer_type> *all_involved_extensions);
    bool explore_branching(kmer_type start_kmer, int start_strand, char *consensus, int &consensus_length, kmer_type previous_kmer);
    bool explore_branching(kmer_type start_kmer, int start_strand, char *consensus, int &consensus_length, kmer_type previous_kmer, set<kmer_type> *all_involved_extensions);
    bool validate_consensuses(set<string> consensuses, char *result, int &result_length);

public:
    static bool all_consensuses_almost_identical(set<string> consensuses);
    MonumentTraversal(Bloom *given_bloom, Set *given_debloom, Terminator *given_terminator) : Traversal(given_bloom,given_debloom,given_terminator) {}

	bool find_starting_kmer(kmer_type kmer, kmer_type &starting_kmer);
    
};

#endif
