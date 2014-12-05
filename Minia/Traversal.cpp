#include <stdio.h>
#include <algorithm> // for max()
#include "Traversal.h"

using namespace std;

Traversal::~Traversal()
{
}


void Traversal::set_maxlen(int given_maxlen)
{
    maxlen = given_maxlen;
}

void Traversal::commit_stats()
{
    final_stats = stats;
}

void Traversal::revert_stats()
{
    stats = final_stats;
}


// --------------------------
// generic structure for traversals

void Traversal::set_max_depth(int given_max_depth)
{
    max_depth = given_max_depth;
}

void Traversal::set_max_breadth(int given_max_breadth)
{
    max_breadth = given_max_breadth;
}

// mark recorded extensions
void Traversal::mark_extensions(set<kmer_type> *extensions_to_mark)
{
    if (terminator == NULL)
        return;
    for(set<kmer_type>::iterator it = extensions_to_mark->begin(); it != extensions_to_mark->end() ; ++it)
    {
        terminator->mark(*it); 
    }
}

// return the number of extension possibilities, the variable nt will contain one of them 
// order=0: just examine immediate neighbors
// order>0: don't return extensions yielding to deadends of length <= order
// todo-order>0: there is probably a minor bug: it is likely to return 0 extensions too early for end of contigs 
int traversal_extensions(kmer_type kmer, int strand, int &nt, Bloom *bloom_solid_kmers, Set *debloom)
{
    if (order==0) // faster order=0 extensions (note: in fact, order is always equal to 0)
    {
        int nb_extensions = 0;
        for(int test_nt=0; test_nt<4; test_nt++) 
        {
            int current_strand = strand;
            kmer_type current_kmer = next_kmer(kmer,test_nt,&current_strand);

            if (bloom_solid_kmers->contains(current_kmer) && !debloom->contains(current_kmer))
            {
                nt = test_nt;
                nb_extensions ++;
            }
        }
        return nb_extensions;
    }
    else
    {
        if (order==1) // optimized code for order=1 (copied from assemb.cpp)
        {
            int nb_extensions = 0;
            for(int test_nt=0; test_nt<4; test_nt++) 
            {
                int current_strand = strand;
                kmer_type current_kmer = next_kmer(kmer,test_nt, &current_strand);

                if(bloom_solid_kmers->contains(current_kmer) && !debloom->contains(current_kmer)){

                    bool is_linked = false;
                    for(int tip_nt=0; tip_nt<4; tip_nt++) 
                    {
                        int new_strand = current_strand;
                        kmer_type kmer_after_possible_tip = next_kmer(current_kmer,tip_nt, &new_strand);
                        if(bloom_solid_kmers->contains(kmer_after_possible_tip) && !debloom->contains(kmer_after_possible_tip))
                        {
                            is_linked = true;
                            break;
                        }
                    }
                    if (!is_linked)
                        continue; // it's a tip, because it's linked to nothing

                    nt = test_nt;
                    nb_extensions++;
                }
            }
            return nb_extensions;
        }
        else
        { // slower, general code for order>=0
            Frontline frontline( kmer, strand, bloom_solid_kmers, debloom, NULL, 0);
            while (frontline.depth <= order) // go one step further than order
            {
                frontline.go_next_depth();
                if (frontline.size() <= 1) // stop when no more ambiguous choices
                    break;
                if (frontline.size() > 10) // don't allow a breadth too large anyway
                    break;
            }
            if (frontline.size() > 0) // recover the nt that lead to this node
                nt = frontline.front().nt;
            return frontline.size();
        }
    }
}

int Traversal::extensions(kmer_type kmer, int strand, int &nt)
{
    return traversal_extensions(kmer,strand,nt,bloom,debloom);
}

/*
 * this function is actually only used when order > 0
 * tip is: 
 * one strand: 0 extension
 * other strand: 1 extension to N
 * N: >= 3 neighbors
 * */
bool is_tip(kmer_type kmer, Bloom *bloom_solid_kmers, Set *debloom)
{
    int nb_extensions[2]={0};
    kmer_type N=0;
    int N_strand=0;
    for (int strand=0;strand<2; strand++)
    {
        for(int test_nt=0; test_nt<4; test_nt++) 
        {
            int current_strand = strand;
            kmer_type current_kmer = next_kmer(kmer,test_nt, &current_strand);
            if(bloom_solid_kmers->contains(current_kmer) && !debloom->contains(current_kmer)){
                N = current_kmer;
                N_strand = current_strand;
                nb_extensions[strand]++;
            }
        }
    }
    /*    if (nb_extensions[0] != 0 &&  nb_extensions[1] != 0)
          return false;*/
    // fixme-order>0: too strict, because includes ends of contigs
    if (nb_extensions[0] == 0 ||  nb_extensions[1] == 0)
        return true;

    // now test degree of N
    int N_degree = 0;
    for(int test_nt=0; test_nt<4; test_nt++) 
    {
        int current_strand = N_strand;
        kmer_type current_kmer = next_kmer(N,test_nt, &current_strand);

        if(bloom_solid_kmers->contains(current_kmer) && !debloom->contains(current_kmer)){
            N_degree++;
        }
    } 
    return N_degree>=3;

}

// from a branching kmer, get a new node that has never been used before
// (very simple initial k-mer selection, current used in minia-graph)
bool Traversal::get_new_starting_node(kmer_type branching_kmer, kmer_type &starting_kmer)
{
    char newNT[2];
    int nt;

    // start with the branching kmer itself
    if ( ! terminator->is_marked(branching_kmer) )
    {
        terminator->mark(branching_kmer);
        starting_kmer = branching_kmer; 
        return true;
    }

    for (int current_strand = 0; current_strand<2 ; current_strand++)
    {
        for(nt=0; nt<4; nt++) 
        {
            int strand = current_strand;
            kmer_type current_kmer;
            current_kmer = next_kmer(branching_kmer,nt,&strand);

            if (bloom->contains(current_kmer) && !debloom->contains(current_kmer))
            {
                // only start from an unmarked nt/strand combo
                if (terminator->is_marked(current_kmer))
                    continue;

                terminator->mark(current_kmer);

                starting_kmer = current_kmer; 
                return true;
            }
        }
    }

    return false;
}

// improved version of the code above
// TODO: port minia-graph to use this version, make sure it doesn't break the simple paths traversal
bool Traversal::get_new_starting_node_improved(kmer_type branching_kmer, kmer_type &starting_kmer)
{
    char newNT[2];
    int nt;

    for (int current_strand = 0; current_strand<2 ; current_strand++)
    {
        for(nt=0; nt<4; nt++) 
        {
            int strand = current_strand;
            kmer_type current_kmer;
            current_kmer = next_kmer(branching_kmer,nt,&strand);

            if (bloom->contains(current_kmer) && !debloom->contains(current_kmer))
            {
                // alright let's use this convention now:
                // to select new kmers: mark non-branching neighbors of branching nodes
                // to mark as used in assembly: mark branching nodes

                // only start from a non-branching k-mer
                if (terminator->is_branching(current_kmer))
                    continue;

                if (terminator->is_marked(current_kmer))
                    continue;
                terminator->mark(current_kmer);

                starting_kmer = current_kmer; 
                return true;
            }
        }
    }

    // actually never start in a branching k-mer
//    if ( ! terminator->is_marked(branching_kmer) )
//    {
//        terminator->mark(branching_kmer);
//        starting_kmer = branching_kmer; 
//        return true;
//    }
 
    return false;
}



/* the initial k-mer selection function from original Minia release (up until fall 2013) was:
        detect a 2k+2 simple path (anything NOT deadend or snp) around the branching kmer and start to extend from it
   the rationale was:
         -> what if it's a simplepathtraversal and it chooses to start in a deadend?
         -> what if it's a monumenttraversal and the kmer is a true branching: won't be traversed
         yet, branching kmers are the only indexed ones
 //bool Traversal::find_starting_kmer_inside_simple_path(kmer_type branching_kmer, kmer_type &starting_kmer)
 * this function was unused since fall 2013 and to prevent bugs due to new marking conventions, it has been removed on 10 Feb 2014.
 * now we use this:
  -> get a better starting point than a branching kmer inside a long (~2k) simple path:  
  -> a k-mer that isn't inside a bubble/tip
*/
bool MonumentTraversal::find_starting_kmer(kmer_type branching_kmer, kmer_type &starting_kmer)
{
    char newNT[2];
    int nt;
    bool debug=false;
    int sum_depths = 0;

    if (!get_new_starting_node_improved(branching_kmer, starting_kmer))
        return false;

    if (debug) printf("getting new starting kmer\n");
   
    for (int strand = 0; strand<2 ; strand++)
    {
        kmer_type previous_kmer = 0;
        int previous_strand = 0;

        // do a BFS to make sure we're not inside a bubble or tip
        Frontline frontline( starting_kmer, strand, bloom, debloom, terminator, NULL, 0, false);

        do
        {
            bool should_continue = frontline.go_next_depth();
            if (!should_continue)
            {
                if (debug) printf("strand %d shouldnt continue\n", strand);
                break;
            }

            // put the same contraints as in a bubble
            if (frontline.depth > max_depth || frontline.size() > max_breadth)
            {
                if (debug) printf("strand %d reached max depth or breadth (%d %d)\n",strand, frontline.depth,frontline.size());
                break;
            }

            // stopping condition: nothing more to explore
            if (frontline.size() == 0)
            {
                if (debug) printf("strand %d nothing more to explore\n", strand);
                break;
            }
           
            char useless_string[max_depth+1];
            int useless_int;
            
            if (frontline.size() <= 1)
            {
                kmer_type current_kmer = 0;
                if (frontline.size() == 1)
                {
                    node current_node = frontline.front();
                    current_kmer = current_node.kmer;
                }

                if ((previous_kmer != 0) && terminator->is_branching(previous_kmer))
                {
                    /* the current situation is:
                     *      
                     *    current_kmer  previous_kmer
                     *   -O-------------O------------- ... ---starting_kmer
                     *                  \_....
                     *
                     * or 
                     *
                     *   [no extension] previous_kmer
                         *   X              O------------- ... ---starting_kmer
                     *                  \_....
                     * 
                     *   so, by looking one k-mer ahead, we make sure that previous_kmer only branches to the right
                     *
                     */
                    set<kmer_type> all_involved_extensions;
                    Terminator *save_terminator = terminator;
                    terminator = NULL; // do not use terminator in the following bubble traversal
                    if (explore_branching(previous_kmer, 1-previous_strand, (char*)useless_string, useless_int, current_kmer, &all_involved_extensions))
                    {
                        if (debug) printf("depth %d useless int %d and starting belongs %d nb involved nodes %lu\n",frontline.depth,useless_int,all_involved_extensions.find(starting_kmer) != all_involved_extensions.end(),all_involved_extensions.size());
                        if (all_involved_extensions.find(starting_kmer) != all_involved_extensions.end())
                        {
                            terminator = save_terminator;
                            return false; // starting_kmer is in a tip/bubble starting from current_kmer
                        }
                        
                    }
                    terminator = save_terminator;
                }
            }
            // update previous_kmer
            if (frontline.size() == 1)
            {
                node current_node = frontline.front();
                kmer_type current_kmer = current_node.kmer;
                previous_kmer = current_kmer;
                previous_strand = current_node.strand;
            }
            else
                previous_kmer = 0;
        }
        while (1);

        if (debug) printf("strand %d depth %d\n",strand,frontline.depth);
        sum_depths += frontline.depth;
    }
    
    // don't even assemble those regions which have no chance to yield a long contig
    if (sum_depths < (sizeKmer+1))
        return false;

    return true;
}

// main traversal function, which calls avance()
// important:
// for MonumentTraversal, either "previous_kmer" is unspecified and then "starting_kmer" is required to be non-branching,
// or, if starting_kmer is branching, please specify the "previous_kmer" parameter, corresponding to a left k-mer that will 
// be ignored during in-branching checks
int Traversal::traverse(kmer_type starting_kmer, char* resulting_sequence, int starting_strand, kmer_type previous_kmer)
{
    kmer_type current_kmer = starting_kmer;
    int current_strand = starting_strand; // 0 = forward, 1 = reverse; 
    int len_extension = 0;
    char newNT[max_depth+1];
    int nnt = 0;
    bool looping = false;

    int bubble_start, bubble_end;
    bubbles_positions.clear();

    //printf(" traversing %llX strand:%d\n",starting_kmer,current_strand);

    while( (nnt=avance(current_kmer, current_strand, len_extension == 0, newNT, previous_kmer)))
    {
        if (nnt < 0) // found branching or marked kmers
            break;

        if (nnt > 1) // it's a bubble for sure
            bubble_start = len_extension;

        // keep re-walking the nucleotides we just discovered, to append to consensus and mark kmers as we go
        for (int cur_nt = 0; cur_nt < nnt; cur_nt++)
        {
            resulting_sequence[len_extension]=newNT[cur_nt];
            len_extension++;
            previous_kmer = current_kmer;

            current_kmer = next_kmer(current_kmer,NT2int(newNT[cur_nt]),&current_strand);
#ifndef DONTMARK
            terminator->mark(current_kmer); // mark kmer as used in the assembly 
#endif

            if (current_kmer == starting_kmer) // perfectly circular regions with no large branching can happen (rarely)
                looping = true;
        }

        if (nnt > 1)
        {
            bubble_end = len_extension;
            bubbles_positions.push_back(std::make_pair(bubble_start,bubble_end));
        }

        if (looping)
            break;

        if (len_extension > maxlen)
        {
            //    fprintf(stderr,"max contig len reached \n");
            break;
        }
    }
    resulting_sequence[len_extension]='\0';
    return len_extension;
}


// ----------------
// random branching traversal

char Traversal::random_unmarked_avance(kmer_type kmer, int current_strand, bool first_extension, char * newNT)
{
    char bin2NT[4] = {'A','C','T','G'};

    int nt;
    for(nt=0; nt<4; nt++) //takes first branch we find
    {
        int strand = current_strand;

        kmer_type new_graine = next_kmer(kmer,nt,&strand);

        if(bloom->contains(new_graine) && !debloom->contains(new_graine) && !terminator->is_marked(kmer, nt, current_strand)){
            *newNT = bin2NT[nt];
            return 1;
        }
    }
    return 0;
}

char RandomBranchingTraversal::avance(kmer_type kmer, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer)
{
    return random_unmarked_avance(kmer,current_strand,first_extension,newNT);
}

// ----------------
// simple paths traversal
// invariant: the input kmer has no in-branching.
// returns:
// 1 if a good extension is found
// 0 if a deadend was reached
// -1 if out-branching was detected
// -2 if no out-branching but next kmer has in-branching
int Traversal::simple_paths_avance(kmer_type kmer, int strand, bool first_extension, char * newNT)
{
    char bin2NT[4] = {'A','C','T','G'};

    int nb_extensions = 0, in_branching_degree = 0;
    int good_nt;

    // return the number of possible forward extensions
    nb_extensions = extensions(kmer, strand, good_nt);

    if (nb_extensions == 1)
    {
        // if the next kmer has in-branching, don't extend the current kmer
        int second_strand = strand;
        kmer_type second_kmer = next_kmer(kmer,good_nt,&second_strand);
        int osef;
        in_branching_degree = extensions(second_kmer,1-second_strand,osef);
        if (in_branching_degree > 1) 
            return -2;

        *newNT = bin2NT[good_nt];
        return 1;
    }

    if (nb_extensions > 1) // if this kmer has out-branching, don't extend it
        return -1;

    return 0;
}

char SimplePathsTraversal::avance(kmer_type kmer, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer)
{
    return max(simple_paths_avance(kmer,current_strand,first_extension,newNT),0);
}


// ----------------
// 
// Monument traversal


// a frontline is a set of nodes having equal depth in the BFS
Frontline::Frontline(kmer_type starting_kmer, int starting_strand, Bloom *bloom, Set *debloom, Terminator *terminator, set<kmer_type> *all_involved_extensions, kmer_type previous_kmer, bool check_in_branching) :
    starting_kmer(starting_kmer), starting_strand(starting_strand), bloom(bloom), debloom(debloom), terminator(terminator), all_involved_extensions(all_involved_extensions), previous_kmer(previous_kmer), check_in_branching(check_in_branching), depth(0)
{
    already_frontlined.insert(starting_kmer); 
    already_frontlined.insert(previous_kmer);

    node first_node(starting_kmer,starting_strand,-1);
    frontline.push(first_node);
}

bool Frontline::go_next_depth() 
{
    // extend all nodes in this frontline simultaneously, creating a new frontline
    stopped_reason=NONE;
    queue_nodes new_frontline;
    while (!frontline.empty())
    {
        node current_node = frontline.front();
        frontline.pop();
        kmer_type current_kmer = current_node.kmer;
        int current_strand = current_node.strand;

        // make sure this node doesn't have large in-branching.
        if ( check_in_branching && check_inbranching(current_kmer,current_strand))
        {
            //printf("######## found large in-branching (depth=%d)\n",depth);
            return false; // detected that the bubble isn't simple (there is large in-branching inside)
        }

        // enqueue all neighbors of this node, except these that were already in a frontline
        for(int nt=0; nt<4; nt++)
        {
            kmer_type new_kmer = current_kmer;
            int new_strand = current_strand; 

            // propagate information where this node comes from
            int from_nt = (current_node.nt == -1) ? nt : current_node.nt;

            // go to next kmer            
            new_kmer = next_kmer(new_kmer,nt,&new_strand);
           
            // test if that node hasn't already been explored
            if (already_frontlined.find(new_kmer)!= already_frontlined.end())
                continue;

            if(bloom->contains(new_kmer) && ((debloom == NULL) || (!debloom->contains(new_kmer))))
            {

                // if this bubble contains a marked (branching) kmer, stop everyone at once (to avoid redundancy)
                if (terminator != NULL && terminator->is_branching(new_kmer))
                    if (terminator->is_marked_branching(new_kmer))
                    {
                        stopped_reason=Frontline::MARKED;
                        return false; 
                    }

                node new_node(new_kmer,new_strand,from_nt);
                new_frontline.push(new_node);
                already_frontlined.insert(new_kmer);

                //if (check_in_branching)
                    //printf("frontline _depth: %d enqueuing kmer %s\n",depth ,print_kmer(new_kmer));

                // since this extension is validated, insert into the list of involved ones 
                if (all_involved_extensions != NULL)
                    all_involved_extensions->insert(new_kmer);
            }
        }
    }
    frontline = new_frontline;
    ++depth;
    return true;
}

int Frontline::size()
{
    return frontline.size();
}

node Frontline::front()
{
    return frontline.front();
}

// new code, not in monument, to detect any in-branching longer than 3k
bool Frontline::check_inbranching(kmer_type from_kmer, int from_strand)
{
    int nt;

    for(nt=0; nt<4; nt++) 
    {
        int strand = 1-from_strand;
        kmer_type current_kmer;
        current_kmer = next_kmer(from_kmer,nt,&strand);

        // only check in-branching from kmers not already frontlined 
        // which, for the first extension, includes the previously traversed kmer (previous_kmer)
        // btw due to avance() invariant, previous_kmer is always within a simple path
        if (already_frontlined.find(current_kmer) != already_frontlined.end())
            continue;

        if (bloom->contains(current_kmer) && !debloom->contains(current_kmer))
        {

            // create a new frontline inside this frontline to check for large in-branching (i know, we need to go deeper, etc..)
            Frontline frontline( current_kmer, strand, bloom, debloom, terminator, all_involved_extensions,from_kmer, false);

            do
            {
                bool should_continue = frontline.go_next_depth();
                if (!should_continue)
                {
                    stopped_reason=Frontline::IN_BRANCHING_OTHER;
                    break;
                }
                // don't allow a depth > 3k
                if (frontline.depth > 3 * sizeKmer)
                {
                    stopped_reason=Frontline::IN_BRANCHING_DEPTH;
                    break;
                }

                // don't allow a breadth too large
                if (frontline.size()> 10)
                {
                    stopped_reason=Frontline::IN_BRANCHING_BREADTH;
                    break;
                }

                // stopping condition: no more in-branching
                if (frontline.size() == 0)
                    break;
            }
            while (1);

            if (frontline.size() > 0)
                return true; // found large in-branching
        }
    }

    // didn't find any in-branching
    return false;
}

// similar to Monument's extension_graph.py:find_end_of_branching
// basically do a bounded-depth, bounded-breadth BFS
int MonumentTraversal::find_end_of_branching(kmer_type starting_kmer, int starting_strand, kmer_type &end_kmer, int &end_strand, kmer_type previous_kmer, set<kmer_type> *all_involved_extensions)
{
    bool check_in_branching = true;
    Frontline frontline( starting_kmer, starting_strand, bloom, debloom, terminator, all_involved_extensions, previous_kmer, check_in_branching);

    do
    {
        bool should_continue = frontline.go_next_depth();
        if (!should_continue)
        {
            if (frontline.stopped_reason == Frontline::MARKED)
                stats.couldnt_because_marked_kmer++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_DEPTH)
                stats.couldnt_inbranching_depth++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_BREADTH)
                stats.couldnt_inbranching_breadth++;
            if (frontline.stopped_reason == Frontline::IN_BRANCHING_OTHER)
                stats.couldnt_inbranching_other++;
            return 0;
        }

        // don't allow a depth too large
        if (frontline.depth > max_depth)
        {
            stats.couldnt_traverse_bubble_depth++;
            return 0;
        }

        // don't allow a breadth too large
        if (frontline.size()> max_breadth)
        {
            stats.couldnt_traverse_bubble_breadth++;
            return 0;
        }

        // stopping condition: frontline is either empty, or contains only 1 kmer
        // needs the kmer to be non-branching, in order to avoid a special case of bubble immediatly after a bubble
        // affects mismatch rate in ecoli greatly
        if (frontline.size() == 0)
        {
            stats.couldnt_find_extension++;
            return 0;
        }
        
//      if (frontline.size() == 1) // longer contigs but for some reason, higher mismatch rate
        if (frontline.size() == 1 && (terminator == NULL || ( !terminator->is_branching(frontline.front().kmer) )))
            break;
    }
    while (1);

    if (frontline.size()==1)
    {
        node end_node = frontline.front();
        end_kmer = end_node.kmer;
        end_strand = end_node.strand;
        return frontline.depth;
    }

    return 0;
}

// similar to Monument's extension_graph.py:all_paths_between
set<string> MonumentTraversal::all_consensuses_between(kmer_type start_kmer, int start_strand, kmer_type end_kmer, int end_strand, int traversal_depth, set<kmer_type> used_kmers, string current_consensus, bool &success)
{
    char bin2NT[4] = {'A','C','T','G'};
    //printf("all consensuses between traversal_depth: %d kmer %s success %d\n",traversal_depth,print_kmer(start_kmer),success);
    set<string> consensuses;

    // find_end_of_branching and all_consensues_between do not always agree on clean bubbles ends
    // until I can fix the problem, here is a fix
    // to reproduce the problem: SRR001665.fasta 21 4
    if (traversal_depth < -1)
    {
        success = false;
        return consensuses; 
    }

    if (start_kmer == end_kmer)// not testing for end_strand anymore because find_end_of_branching doesn't care about strands
    {
        consensuses.insert(current_consensus);
        return consensuses;
    }

    // visit all neighbors
    for(int nt=0; nt<4; nt++)
    {
        // craft neighbor node
        int new_strand = start_strand;
        kmer_type new_graine = next_kmer(start_kmer,nt,&new_strand);

        // check if the neighbor node is valid 
        if(bloom->contains(new_graine) && !debloom->contains(new_graine)){

            // don't resolve bubbles containing loops 
            // (tandem repeats make things more complicated)
            // that's a job for a gapfiller
            if (used_kmers.find(new_graine) != used_kmers.end())
            {
                success = false;
                return consensuses;
            }

            // generate extended consensus sequence
            string extended_consensus(current_consensus);
            extended_consensus.append(1,bin2NT[nt]);

            // generate list of used kmers (to prevent loops)
            set<kmer_type> extended_kmers(used_kmers);
            extended_kmers.insert(new_graine);

            // recursive call to all_consensuses_between
            set<string> new_consensuses = all_consensuses_between(new_graine, new_strand, end_kmer, end_strand, traversal_depth - 1, extended_kmers, extended_consensus, success);
            consensuses.insert(new_consensuses.begin(), new_consensuses.end());

            // mark to stop we end up with too many consensuses
            if (consensuses.size() > (unsigned int )max_breadth)
                success = false;
        }

        // propagate the stop if too many consensuses reached
        if (success == false)
            return consensuses;
    }
    return consensuses;
}

// just a wrapper
set<string> MonumentTraversal::all_consensuses_between(kmer_type start_kmer, int start_strand, kmer_type end_kmer, int end_strand, int traversal_depth, bool &success)
{
    set<kmer_type> used_kmers;
    used_kmers.insert(start_kmer);
    string current_consensus;
    success = true;
    //printf("all cons between - end kmer = %s\n",print_kmer(end_kmer));
    return all_consensuses_between(start_kmer, start_strand, end_kmer, end_strand, traversal_depth, used_kmers, current_consensus, success);
}

// similar to Monument's extension_graph.py:validate_paths
// return true if, basically, the consensuses aren't too different
bool MonumentTraversal::validate_consensuses(set<string> consensuses, char *result, int &result_length)
{
    bool debug = false;
    // compute mean and stdev of consensuses
    int mean = 0;
    int path_number = 0;
    for(set<string>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        if (debug)  printf("bubble path %d: %s (len=%lu)\n",path_number,(*it).c_str(),(*it).length());
        mean+=(*it).length();
        path_number++;
    }
    mean/=consensuses.size();
    double stdev = 0;
    for(set<string>::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        int consensus_length = (*it).length();
        stdev += pow(fabs(consensus_length-mean),2);
    }
    stdev = sqrt(stdev/consensuses.size());

    // don't traverse large bubbles
    if (mean > max_depth) 
        return false;

    // don't traverse large deadends (here, having one consensus means the other paths were large deadends)
    if (consensuses.size() == 1 && mean > sizeKmer+1) // deadend length should be < k+1 (most have length 1, but have seen up to 10 in ecoli)
        return false;

    if (debug) printf("%lu-bubble mean %d, stdev %.1f\n",consensuses.size(),mean,stdev);

    // traverse bubbles if paths have roughly the same length
    if (stdev>mean/5)
        return false;

    // check that all consensuses are similar
    if (! all_consensuses_almost_identical(consensuses))
        return false;

    // if all good, an arbitrary consensus is chosen
    string chosen_consensus = *consensuses.begin();
    result_length = chosen_consensus.length();
    if  (result_length> max_depth) // it can happen that consensus is longer than max_depth, despite that we didn't explore that far (in a messy bubble with branchings inside)
        return false;

    chosen_consensus.copy(result, result_length);
    return true;
}

bool MonumentTraversal::all_consensuses_almost_identical(set<string> consensuses)
{
    for (set<string>::iterator it_a = consensuses.begin(); it_a != consensuses.end(); it_a++)
    {
        set<string>::iterator it_b = it_a;
        advance(it_b,1);
        while (it_b != consensuses.end())
        {
            if (needleman_wunch(*it_a,*it_b) * 100 < consensuses_identity)
                return false;
            advance(it_b,1);
        }
    }
    return true;
}


// similar to Monument's extension_graph.py:explore_branching
// return true if the branching can be traversed, and mark all involved nodes
bool MonumentTraversal::explore_branching(kmer_type start_kmer, int start_strand, char *consensus, int &consensus_length, kmer_type previous_kmer, set<kmer_type> *all_involved_extensions)
{
    kmer_type end_kmer;
    int end_strand;

    // find end of branching, record all involved extensions (for future marking)
    // it returns false iff it's a complex bubble
    int traversal_depth = find_end_of_branching(start_kmer, start_strand, end_kmer, end_strand, previous_kmer, all_involved_extensions);
    if (!traversal_depth)
    {   
        stats.couldnt_find_all_consensuses++;
        return false;
    }

    // find all consensuses between start node and end node 
    set<string> consensuses;
    bool success;
    consensuses = all_consensuses_between(start_kmer, start_strand, end_kmer, end_strand, traversal_depth+1, success);
    // if consensus phase failed, stop
    if (!success)
        return false;

    // validate paths, based on identity
    bool validated = validate_consensuses(consensuses, consensus, consensus_length);
    if (!validated)
    {
        stats.couldnt_validate_consensuses++;
        return false;
    }

    // the consensuses agree, mark all the involved extensions 
    // (corresponding to alternative paths we will never traverse again)
    mark_extensions(all_involved_extensions);

    return true;
}

// wrapper
bool MonumentTraversal::explore_branching(kmer_type start_kmer, int start_strand, char *consensus, int &consensus_length, kmer_type previous_kmer)
{
    set<kmer_type> *all_involved_extensions = new set<kmer_type>;
    bool res = explore_branching(start_kmer, start_strand, consensus, consensus_length, previous_kmer, all_involved_extensions);
    delete all_involved_extensions;
    return res;
}


// invariant here:
// kmer is always obtained after traversing a non-branching kmer
// in other words, the only in-branching of that kmer is previous_kmer
char MonumentTraversal::avance(kmer_type kmer, int current_strand, bool first_extension, char * newNT, kmer_type previous_kmer)
{

    // if we're on a simple path, just traverse it
    int is_simple_path = simple_paths_avance(kmer, current_strand, first_extension, newNT);
    if (is_simple_path > 0)
        return 1;

    // the following function does:
    // * a bfs from the starting kmer, stopping when:
    //      - breadth > max_breadth
    //      or
    //      - depth > max_depth
    // * check if there a single end point
    // * computing all possible paths between start and end
    // * returns one flattened consensus sequence
    int newNT_length;
    bool success = explore_branching(kmer, current_strand, newNT, newNT_length, previous_kmer);
    if (!success)
    {
        stats.ended_traversals++;
        return 0;
    }

    return newNT_length;
}
