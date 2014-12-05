#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <chrono>

#include "lm.h"
#include "input.h"

/*
 * Constuct the compacted de bruijn graph from list of distinct kmers
 */

using namespace std;

// 10 chars to store the integer hash of a minimizer
#define MINIMIZER_STR_SIZE 10

// keep largest bucket seen, disregarding small ones
unsigned long nb_elts_in_largest_bucket = 10000;

string minimizer2string(int input_int)
{
    long long i = input_int;
    string str = to_string(i);
    assert(str.size() <= MINIMIZER_STR_SIZE);
    if(MINIMIZER_STR_SIZE > str.size())
        str.insert(0, MINIMIZER_STR_SIZE- str.size(), '0');
    return str;
}

bool nextchar(char * c){
	switch(*c){
		case 'a':
		*c='c';
		return true;
		case 'c':
		*c='g';
		return true;
		case 'g':
		*c='t';
		return true;
		case 't':
		*c='a';
		return false;
		default :
		cout<<"Problem with nextchar: "<<c;
		assert(0);
		return false;
	}
}


// reads k-mers and Put kmers in superbuckets
// (other behavior: just count m-mers)
void sortentry(string namefile, const int k, const int m, bool create_buckets, bool m_mer_count){
	int numbersuperbucket(pow(4,m));
	ofstream out[2100];
    ofstream checkpoint;
    InputDot ind;
    if (create_buckets)
    {
    	for(long long i(0);i<numbersuperbucket;i++)
    		out[i].open(".bcalmtmp/z"+to_string(i),ofstream::app);
    }

    ind.init_input(namefile,k);

    while (1)
    {
        string kmer = ind.next_input(k);

        if (kmer == "")
            break;

        if (m_mer_count)
        {
            count_m_mers(kmer, 2*m, k);
            continue;
        }

        int middlemin, leftmostmin, rightmostmin, min;
        
        middlemin=minimiserrc(kmer.substr(1,k-2),2*m);
        leftmostmin=minimiserrc(kmer.substr(0,2*m),2*m);
        rightmostmin=minimiserrc(kmer.substr(kmer.size()-2*m,2*m),2*m);

        int leftmin = (leftmostmin < middlemin) ? leftmostmin : middlemin;

        int rightmin = (rightmostmin < middlemin) ? rightmostmin : middlemin ;

        min = (leftmin < rightmin) ? leftmin : rightmin;

        uint64_t h = min/numbersuperbucket;

        if (create_buckets)
            out[h]<<kmer<< minimizer2string(leftmin) << minimizer2string(rightmin) <<";";
        else
            cout << h << ":" << kmer<< minimizer2string(leftmin) << minimizer2string(rightmin) << ";\n";

        //checkpoint << kmer << minimizer2string(leftmin) << minimizer2string(rightmin) <<";\n";
}
    if (create_buckets)
        cout << "initial partitioning done" << endl;
}



//copy n characters from in to out
void copylm(ifstream* in,int64_t n,ofstream* out){
	int buffsize(1000000),nbbuffer(n/buffsize);
	string str;
	for(int j(0);j<nbbuffer;j++)	{
		str=readn(in,buffsize);
		*out<<str;
	}
	str=readn(in,n-nbbuffer*buffsize);
	*out<<str;
}

void copylmrv(ifstream* in,int64_t n,ofstream* out){
	int buffsize(1000000),nbbuffer(n/buffsize);
	int64_t pos(in->tellg());
	string str;
	int j;
	for(j=1;j<=nbbuffer;j++){
		in->seekg(pos+n-j*buffsize,ios::beg);
		str=readn(in,buffsize);
		*out<<reversecompletment(str);
	}
	int64_t rest=n-nbbuffer*buffsize;
	if(rest!=0){
		in->seekg(pos,ios::beg);
		str=readn(in,rest);
		*out<<reversecompletment(str);
	}
}


//Put nodes from superbuckets to buckets
void createbucket(const string superbucketname,const int m){
	int superbucketnum(stoi(superbucketname));
	ifstream in(".bcalmtmp/z"+superbucketname);
	if(!in.is_open()){
		cerr<<"Problem with Createbucket"<<endl;
		return;
	}
	in.seekg(0, ios_base::end);
	int64_t size(in.tellg()), buffsize(1000000), numberbuffer(size/buffsize),nb(pow(4,m)),suffix;
	if(size==0){
		remove((".bcalmtmp/z"+superbucketname).c_str());
		return;
	}
	in.seekg(0,ios::beg);
	ofstream out[5000];
	for(long long i(0);i<nb;i++){
		out[i].open(".bcalmtmp/"+to_string(superbucketnum*nb+i),ofstream::app);
	}
	int64_t lastposition(-1),position(0),point(0),mini;
	string buffer;
	vector<string> miniv;
	in.seekg(0);

	for(int j(0);j<=numberbuffer;j++){
		if(j==numberbuffer)
        {
			if(size-numberbuffer*buffsize-1!=-1){
				buffer=readn(&in,size-numberbuffer*buffsize-1);
				buffer+=";";
			}
			else
				buffer="";
        }
		else
		{
			buffer=readn(&in,buffsize);
			point+=buffsize;
		}
		for (uint64_t i(0); i<buffer.size(); i++,position++)
        {
			if((buffer)[i]==';'){
                int leftmin, rightmin;
				if(i>=(uint64_t)(MINIMIZER_STR_SIZE * 2))
                {
					leftmin = stoi(buffer.substr(i-(MINIMIZER_STR_SIZE * 2),MINIMIZER_STR_SIZE));
					rightmin = stoi(buffer.substr(i-MINIMIZER_STR_SIZE,MINIMIZER_STR_SIZE));
                }
				else{
					in.seekg(position-(MINIMIZER_STR_SIZE * 2));
					leftmin = stoi(readn(&in,MINIMIZER_STR_SIZE));
					rightmin = stoi(readn(&in,MINIMIZER_STR_SIZE));
                }
				string first_bucket_of_superbucket(to_string((long long)(superbucketnum*nb-1)));
				mini = minbutbiggerthan(leftmin, rightmin, first_bucket_of_superbucket);
				suffix=mini%nb;
				in.seekg(lastposition+1,ios_base::beg);
				copylm(&in,position-lastposition,&out[suffix]);
				lastposition=position;
			}
        }
		in.seekg(point);
	}
	remove((".bcalmtmp/z"+superbucketname).c_str());
}



//count the length of each node
vector<int64_t> countbucket(const string& name){
	vector<int64_t> count;
	ifstream in(".bcalmtmp/"+name);
	if(in){
		in.seekg( 0 , ios_base::end );
		int64_t size(in.tellg()),buffsize(10),numberbuffer(size/buffsize),lastposition(-1),position(0);
		if(size<2)
			return count;
		in.seekg(0,ios::beg);
		string buffer;

		for(int j(0);j<numberbuffer;j++){
			buffer=readn(&in,buffsize);
			for (uint64_t i(0); i<buffer.size(); i++,position++)
				if((buffer)[i]==';'){
					count.push_back(position-lastposition-1);
					lastposition=position;
				}
		}

		buffer=readn(&in,size-numberbuffer*buffsize);
		for (uint64_t i(0); i<buffer.size(); i++,position++)
			if((buffer)[i]==';'){
				count.push_back(position-lastposition-1);
				lastposition=position;
			}
	}
	return count;
}



//true iff node does not contain tag
bool notag(const string& node,const int64_t start,int64_t* n){
	for(uint64_t i(start);i<node.size();i++){
		if (node[i] >= '0' && node[i] <= '9'){
			*n=i;
			return false;
		}
	}
	return true;
}



//return length of tag
int taglength(const string& node, int64_t j){
	int n=1;
	for(uint64_t i(j+1);i<node.size();i++)
		if ((node[i] >= '0' && node[i] <= '9') || node[i]=='+' || node[i]=='-')
			n++;
		else
			return n;
	return n;
}



//Write a node remplacing tags by their sequences
void writeit(const string& outfile,const string& node, int leftmin, int rightmin, vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,int64_t j,const string& fout){
	ofstream out(".bcalmtmp/"+outfile,ios::app);
	char rc;
	if(out){
		int64_t lastposition(0),tag,tagl,position,length;
		pair<int64_t,int64_t> pair;
		do{
			out<<node.substr(lastposition,j-lastposition-1);
			tagl=taglength(node,j);
			rc=node[j-1];
			if(rc=='+'){
				tag=stoi(node.substr(j,tagl));
			}
			else{
				tag=stoi(reversecompletment(node.substr(j,tagl)));
			}
			lastposition=j+tagl;
			pair=(*tagsposition)[tag];
			position=pair.first;
			length= pair.second;
			tagfile->seekg(position,ios_base::beg);
			if(rc=='+')
				copylm(tagfile,length,&out);
			if(rc=='-')
				copylmrv(tagfile,length,&out);
		}
		while(!notag(node,lastposition,&j));
		if(outfile!=fout)
			out<<node.substr(lastposition)<< minimizer2string(leftmin) << minimizer2string(rightmin) <<";";
		else
			out<<node.substr(lastposition) <<";"<<endl;
	}
	else
		cerr<<"writeitbug"<<endl;
}

void put(const string& outfile,const string& node, int leftmin, int rightmin, const string& fout){
	ofstream out(".bcalmtmp/"+outfile,ios::app);
	if(outfile==fout)
	    out<<node <<";" << endl;
    else
	    out<<node<< minimizer2string(leftmin) << minimizer2string(rightmin) <<";";
}

void putorwrite(const string& outfile, const string& node, int leftmin, int rightmin, vector<pair<int64_t,int64_t>>* tagsposition , ifstream* tagfile,const string& fout){
	int64_t i;
	if(notag(node,0,&i))
		put(outfile,node,leftmin, rightmin, fout);
	else
		writeit(outfile,node, leftmin, rightmin, tagsposition,tagfile,i,fout);
}

//Decide where to put a node
void goodplace(const string& node, int leftmin, int rightmin, const string& bucketname,vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,const int m,const string& nameout){
	int nb(pow(4,m)),prefixnumber(stoi(bucketname)/nb+1);
	long long mini(minbutbiggerthan(leftmin, rightmin, bucketname));
	if(mini==-1)
		putorwrite( nameout, node, leftmin, rightmin, tagsposition,tagfile,nameout);
	else{
		long long minipre(mini/nb);
		string miniprefix('z'+to_string(minipre));
		putorwrite( ((minipre >= prefixnumber) ?  miniprefix: to_string(mini)), node, leftmin, rightmin, tagsposition,tagfile,nameout);
	}
}



//Compact a bucket and put the nodes on the right place
void compactbucket(const int& prefix,const int& suffix,const int k,const char *nameout,const int m){
	int64_t buffsize(k),postags(0),length,nb(pow(4,m));
	long long tagnumber(0),numberbucket(prefix*nb+suffix);
	string fullname(to_string(numberbucket)),node,tag,end;
	auto count(countbucket(fullname));
	if(count.size()==0)	{
		remove((".bcalmtmp/"+fullname).c_str());
		return;
	}

    // keep largest bucket seen
    if (count.size() > nb_elts_in_largest_bucket)
    {
        nb_elts_in_largest_bucket = count.size();
	    system(("cp .bcalmtmp/" + fullname + " largest_bucket.dot").c_str());
    }

	ifstream in(".bcalmtmp/"+fullname);
	ofstream tagfile(".bcalmtmp/tags"),out(".bcalmtmp/"+(string)nameout,ios_base::app);
	graph g(k);
	vector<pair<int64_t,int64_t>> tagsposition;

    // add nodes to graph
	if(in && tagfile && out){
		for(auto it=count.begin();it!=count.end();it++){
			length=*it;
			if(length-(MINIMIZER_STR_SIZE*2)<=2*buffsize){
				node=readn(&in,length+1);
				g.addvertex(node.substr(0,length-(MINIMIZER_STR_SIZE*2)));
                g.addleftmin(stoi(node.substr(length-(MINIMIZER_STR_SIZE*2),MINIMIZER_STR_SIZE)));
                g.addrightmin(stoi(node.substr(length-MINIMIZER_STR_SIZE,MINIMIZER_STR_SIZE)));
			}
			else{
				node=readn(&in,buffsize);
				tag=to_string(tagnumber);
				node+="+"+tag+"+";
				tagnumber++;
				copylm(&in,length-(MINIMIZER_STR_SIZE*2)-2*buffsize,&tagfile);
				tagsposition.push_back(make_pair(postags,length-(MINIMIZER_STR_SIZE*2)-2*buffsize));
				postags+=length-(MINIMIZER_STR_SIZE*2)-2*buffsize;
				end=readn(&in,buffsize);
				node+=end.substr(0,buffsize);
				g.addvertex(node);
                g.addleftmin(stoi(readn(&in,MINIMIZER_STR_SIZE)));
                g.addrightmin(stoi(readn(&in,MINIMIZER_STR_SIZE)));
                readn(&in,1); // the ';'
			}
		}
		tagfile.close();
	}

	remove((".bcalmtmp/"+fullname).c_str());
	g.debruijn();

	g.compressh(stoi(fullname));
	ifstream fichiertagin(".bcalmtmp/tags");
    int node_index = 0;
    
	for(auto it(g.nodes.begin());it!=g.nodes.end();it++)
    {
		if(it->size()!=0)
        {
            int leftmin = g.leftmins[node_index];
            int rightmin = g.rightmins[node_index];
			goodplace(*it, leftmin, rightmin,fullname,&tagsposition,&fichiertagin,m,nameout);
        }
        node_index++;
    }

	remove(".bcalmtmp/tags");
	return;
}

//Create a file with the nodes of the compacted graph
void createoutfile(const char *namein,const char *nameout,const int k,const int m){
	auto start=chrono::system_clock::now();
	HashMap hm(build_hash_map(2*m));
	int64_t nbsuperbucket(pow(4,m)),sys(0);
	mkdir(".bcalmtmp",0777);
	sys+=system("rm -rf .bcalmtmp/*");
	remove(nameout);

    // create the hash function
    init_m_mers_table(2*m);
    sortentry(namein,k,m, false, true);
    create_hash_function_from_m_mers(2*m);

	sortentry(namein,k,m);
	for(long long i(0);i<nbsuperbucket;i++){
		createbucket(to_string(i),m);
		for(int j(0);j<pow(4,m);j++)
			compactbucket(i,j,k,nameout,m);
	}
	sys+=system(("mv .bcalmtmp/"+(string)nameout+" "+(string)nameout).c_str());
	if(sys!=0)
		cerr<<"system call failed"<<endl;
	 auto end=chrono::system_clock::now();
	 auto waitedFor=end-start;
	 cout<<"Last for "<<chrono::duration_cast<chrono::seconds>(waitedFor).count()<<" seconds"<<endl;
}

