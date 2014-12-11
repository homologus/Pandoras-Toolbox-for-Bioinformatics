#include "GraphOutput.h"
const string GraphOutput::graph_file_suffix = ".graph";
const string GraphOutput::nodes_file_suffix = ".nodes";
const string GraphOutput::edges_file_suffix = ".edges";
const string GraphOutput::xml_file_suffix = ".xgmml";
const string GraphOutput::json_nodes_file_suffix = ".json_nodes";
const string GraphOutput::json_edges_file_suffix = ".json_edges";
const string GraphOutput::json_file_suffix = ".json";

/************************************************************************************************************************/
/* init function initialize the files need to construct graph file or sequences file					*/
/*															*/
/************************************************************************************************************************/
void GraphOutput::init(bool erase){
  //printf("create a graph erase=%s graph_format=%d  first_id_nodes=%d first_id_edges=%d\n"", erase?"true":"false", graph_format, first_id_els.node, first_id_els.edge);
  graph_file_name=(prefix+graph_file_suffix);
    nodes_file_name=(prefix+nodes_file_suffix);
    edges_file_name=(prefix+edges_file_suffix);
    xml_file_name=(prefix+xml_file_suffix);
    json_nodes_file_name=(prefix+json_nodes_file_suffix);
    json_edges_file_name=(prefix+json_edges_file_suffix);
    json_file_name=(prefix+json_file_suffix);

    switch (graph_format){
        case 0:// FORMAT .GRAPH
            graph_file = fopen(graph_file_name.c_str(),erase?"w":"a");
            fprintf(graph_file,"digraph dedebruijn {\n");
        break;
            
        case 1: // FORMAT .NODES AND .EDGES
            nodes_file = fopen(nodes_file_name.c_str(),erase?"w":"a");
            edges_file = fopen(edges_file_name.c_str(),erase?"w":"a");
        break;
     
        case 2 :// FORMAT .XGMML
            graph_file = fopen(xml_file_name.c_str(),erase?"w":"a");
            //fprintf(graph_file,"<?xml version=\"1.0\"?>\n<!DOCTYPE graph SYSTEM \"http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd\">\n<graph directed=\"1\">\n");
            fprintf(graph_file,"<graph label=\"%s\"\n", (char *)prefix.c_str());
            fprintf(graph_file,"xmlns:dc=\"http://purl.org/dc/elements/1.1/\" \n");
            fprintf(graph_file,"xmlns:xlink=\"http://www.w3.org/1999/xlink\" \n");
            fprintf(graph_file,"xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n");
            fprintf(graph_file,"xmlns:cy=\"http://www.cytoscape.org\" \n");
            fprintf(graph_file,"xmlns=\"http://www.cs.rpi.edu/XGMML\"  \n");
            fprintf(graph_file,"directed=\"1\">\n");
        break;

        case 3: // FORMAT .json
            nodes_file = fopen(json_nodes_file_name.c_str(),erase?"w":"a");
            edges_file = fopen(json_edges_file_name.c_str(),erase?"w":"a");
            graph_file = fopen(json_file_name.c_str(),erase?"w":"a");
        break;
    }
}

/************************************************************************************************************************/
/* 		printf GraphOutput and initialize files (files are not erasing)						*/
/*															*/
/************************************************************************************************************************/
GraphOutput::GraphOutput(string prefix, int graph_format, id_els first_id_els)  : prefix(prefix), graph_format(graph_format), first_id_els(first_id_els)
{ // PIERRE: need something different than 0 for the first node
  printf("graph_format=%d first_id_nodes=%d first_id_edges=%d\n", graph_format, first_id_els.node, first_id_els.edge);
  init(true);
}

/************************************************************************************************************************/
/* 		Initialize first elements and files  (files are erasing)						*/
/*															*/
/************************************************************************************************************************/
GraphOutput::GraphOutput(string prefix, int graph_format) : prefix(prefix), graph_format(graph_format)
{
  first_id_els.node=0;
  first_id_els.edge=0; 
  printf("graph_format=%d first_id_nodes=%d first_id_edges=%d\n", graph_format, first_id_els.node, first_id_els.edge);
  init(true);
}

/************************************************************************************************************************/
/* 		write graph file or sequence file									*/
/*															*/
/************************************************************************************************************************/
void GraphOutput::close()
{
    switch (graph_format){
        
        case 0:
            fprintf(graph_file,"}\n");
            fclose(graph_file);
        break;
    
        case 1:
            fclose(nodes_file);
            fclose(edges_file);
        break;
    
        case 2:
            fprintf(graph_file,"</graph>\n");
            fclose(graph_file);
        break;
      
        case 3:

           // We need to store all nodes and then all edges in the final .json file
            fclose(nodes_file);
            fclose(edges_file);
            ifstream nodes(json_nodes_file_name.c_str(), ios::in);
            ifstream edges(json_edges_file_name.c_str(), ios::in);
            
            if(!edges || !nodes){fprintf(stderr,"Cannot open file %s, %s or %s, exit\n", json_edges_file_suffix.c_str(), json_nodes_file_suffix.c_str()); exit(1);}
            
            string line;

	    fprintf(graph_file,"{\n \"Starter\":[\n{");
	    fprintf(graph_file,"\n \"nodes\": [\n");
	    getline(nodes,line); 
      	    fprintf(graph_file,"%s",line.c_str()); // prints the first node without comma before
	    //for each node
	    while(getline(nodes,line)){
		fprintf(graph_file,",\n%s",line.c_str()); // prints the other nodes
	    };
	    fprintf(graph_file,"\n],\n");
	    fprintf(graph_file,"\"edges\": [\n");
            getline(edges,line);
            fprintf(graph_file,"%s",line.c_str()); // prints the first edge without comma before
	    //for each edge
            while(getline(edges,line)) {
			fprintf(graph_file,",\n%s",line.c_str()); // prints the others edges
            };

	    //end of graph file en close file
	    fprintf(graph_file,"\n]\n}\n");
	    nodes.close(); remove(json_nodes_file_name.c_str());
	    edges.close(); remove(json_edges_file_name.c_str());
            fclose(graph_file);
    }
}


/************************************************************************************************************************/
/* 	recalculate length for a node (more efficient than capture length in string and convert the in integer) 	*/
/*															*/
/************************************************************************************************************************/
long GraphOutput::sequence_length(string line)
{ 
	string  seq_char;
        int err,match,start, end;
        regex_t preg;
        long seq_len=0;
      	size_t nmatch, size;				
	const char *str_regex ="([A-Z]+)"; //regex capture sequences characters
  	const char *line_c =NULL;

	line_c = line.c_str();
  	err = regcomp (&preg, str_regex, REG_EXTENDED);
   				 
	if (err == 0)//security for error string snapshot and if regex match
   	{
		nmatch = 0;
      		nmatch = preg.re_nsub;
		regmatch_t *pmatch=NULL;
     		pmatch = (regmatch_t*) malloc (sizeof (*pmatch) * nmatch);
      		if (pmatch)
      		{
        		match = regexec (&preg, line_c, nmatch, pmatch, 0);
         		regfree (&preg);
         		if (match == 0)
         		{
           			char *seq_char =NULL;
            			start = pmatch[0].rm_so;
            			end = pmatch[0].rm_eo;
            			size = end - start;
                                seq_len = sizeof(line_c[start])*(size);
         		}
         	}
      	}
      	else
      	{
        	fprintf (stderr, "LOW MEMORY !\n");
         	exit (EXIT_FAILURE);
      	} 
return  seq_len;
}

/************************************************************************************************************************/
/* 		output a single node to a file										*/
/*															*/
/************************************************************************************************************************/
void GraphOutput::print_node(long index, char *ascii_node) // output a single node to a file
{
 int len;
    switch (graph_format){
        case 0: // DOT format
            fprintf(graph_file,"%ld [label=\"%s\"];\n",index,ascii_node);
        break;
            
        case 1: // kissplice format
            fprintf(nodes_file,"%ld\t%s\n",index,ascii_node);
        break;
  
        case 2: // XGMML format
            fprintf(graph_file,"<node id=\"%ld\" label=\"%s\">\n</node>\n",index,ascii_node);
        break;

        case 3: // json format
	    string seq = ascii_node;
	    len = seq.size();
            fprintf(nodes_file," { \"data\": { \"id\":\"%ld\", \"length\":%d, \"sequence\":\"%s\"}}\n",index,len,ascii_node);
        break;
  }  
}


/************************************************************************************************************************/
/* 		output a single edges to a file										*/
/*															*/
/************************************************************************************************************************/
void GraphOutput::print_edge(long index, long id, long id2, string label)
{
    switch (graph_format){
        case 0: // DOT format
            fprintf(graph_file,"%ld -> %ld [label=\"%s\"];\n",id,id2,label.c_str());
        break;
  
        case 1: // kissplice format
            fprintf(edges_file,"%ld\t%ld\t%s\n",id,id2,label.c_str());
        break;
  
        case 2: // XGMML format
            fprintf(graph_file,"<edge source=\"%ld\" target=\"%ld\" label=\"%s\">\n</edge>\n",id,id2,label.c_str());
        break;
            
        case 3: // json format
            fprintf(edges_file,"{ \"data\":{ \"id\": \"%ld\", \"source\": \"%ld\",\"target\": \"%ld\",\"direction\": \"%s\"}}\n",index, id,id2,label.c_str());
            //fprintf(edges_file,"{ \"data\":{ \"source\": %ld,\"target\": %ld,\"direction\": \"%s\"}}\n",id,id2,label.c_str());
        break;
  }
  
}

/************************************************************************************************************************/
/* 		load nodes extremities											*/
/*															*/
/************************************************************************************************************************/
void GraphOutput::load_nodes_extremities(string linear_seqs_name)
{
  kmer_links.clear(); // PIERRE: reset previous stored kmer links

    Bank *Nodes = new Bank((char *)linear_seqs_name.c_str());
    long nb_nodes = first_id_els.node; //PIERRE;
    char * rseq;
    int readlen;

    sizeKmer--; // nodes extremities overlap by (k-1)-mers, so let's extract (k-1)-mers

    while (Nodes->get_next_seq(&rseq,&readlen))
    {
        kmer_type left_kmer, right_kmer, left_kmer_fw, left_kmer_rc, right_kmer_fw, right_kmer_rc;
        left_kmer = extractKmerFromRead(rseq,0,&left_kmer_fw,&left_kmer_rc, false);
        right_kmer = extractKmerFromRead(rseq,readlen-sizeKmer,&right_kmer_fw,&right_kmer_rc, false);
        Strand left_strand = (left_kmer == left_kmer_fw)?FW:RC;
        Strand right_strand = (right_kmer == right_kmer_fw)?FW:RC;

        kmer_links[left_kmer].insert(node_strand(nb_nodes, left_strand, LEFT));
        kmer_links[right_kmer].insert(node_strand(nb_nodes, right_strand, RIGHT));

        nb_nodes++;
    }
    Nodes->close();
    delete Nodes;

    sizeKmer++; // make sure to restore k
}


/************************************************************************************************************************/
/* 		construct node file and edge file for graph file							*/
/*															*/
/************************************************************************************************************************/
 id_els GraphOutput::construct_graph(string linear_seqs_name) // PIERRE: i need to know the last nb_nodes
{
    Bank *Nodes = new Bank((char *)linear_seqs_name.c_str());
    id_els nb_els = first_id_els; //Alexan: stucture for print id elements in graph output

    char * rseq;
    int readlen;

    Nodes->rewind_all();

    sizeKmer--; // nodes extremities overlap by (k-1)-mers, so let's extract (k-1)-mers

    // for each node, output all the out-edges (in-edges will correspond to out-edges of neighbors)
    while (Nodes->get_next_seq(&rseq,&readlen))
    {
	
        kmer_type left_kmer, right_kmer, left_kmer_fw, left_kmer_rc, right_kmer_fw, right_kmer_rc;
        set<node_strand>::iterator it;

        left_kmer = extractKmerFromRead(rseq,0,&left_kmer_fw,&left_kmer_rc, false);
        right_kmer = extractKmerFromRead(rseq,readlen-sizeKmer,&right_kmer_fw,&right_kmer_rc, false);
        Strand left_strand = (left_kmer == left_kmer_fw)?FW:RC;
        Strand right_strand = (right_kmer == right_kmer_fw)?FW:RC;


        // left edges (are revcomp extensions)
        for (it = kmer_links[left_kmer].begin(); it != kmer_links[left_kmer].end(); it++)
        {
            long cur_node = it->node;
            Strand cur_strand = it->strand;
            LeftOrRight cur_left_or_right = it->left_or_right;

            if (cur_node ==nb_els.node) // prevent self loops on same kmer
                 if (readlen == sizeKmer)
                    continue;
            
            string label = "R";

            if (cur_left_or_right == LEFT)
            {
                if (cur_strand != left_strand)
                    label+=(string)"F";
                else
                    continue;
            }
            else
            {
                if (cur_strand == left_strand)
                    label+=(string)"R";
                else
                    continue;
            }


            print_edge(nb_els.edge, nb_els.node,cur_node,label);
	        nb_els.edge++; 
        }

        // right edges
        for (it = kmer_links[right_kmer].begin(); it != kmer_links[right_kmer].end(); it++)
        {
            long cur_node = it->node;
            Strand cur_strand = it->strand;
            LeftOrRight cur_left_or_right = it->left_or_right;

            if (cur_node == nb_els.node) // prevent self loops on same kmer
                 if (readlen == sizeKmer)
                    continue;
           
            string label = "F";

            if (cur_left_or_right == LEFT)
            {
                if (cur_strand == right_strand)
                    label+=(string)"F";
                else
                    continue;
            }
            else
            {
                if (cur_strand != right_strand)
                    label+=(string)"R";
                else
                    continue;
            }

            print_edge(nb_els.edge, nb_els.node,cur_node,label);
	        nb_els.edge++;
        }

        //nodes
        print_node(nb_els.node, rseq);   

        nb_els.node++;
    }

    sizeKmer++; // make sure to restore k
    Nodes->close();
    delete Nodes;
    return nb_els;
}


