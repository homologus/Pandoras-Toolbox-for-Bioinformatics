#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include <unordered_map>

using namespace std;

void create_hash_function_from_m_mers(int m);
void count_m_mers(string str, int m, int k);
void init_m_mers_table(int m);

typedef unordered_map<string,int> HashMap;

HashMap build_hash_map(int len);

int shash(const string& s, int& previous_shash, unsigned int start_pos = 0, int length = -1);

string inverse_shash (int num, int len);

int minimiserrc(const string &node,const int &minimisersize);

int minbutbiggerthan(int leftmin, int rightmin, const string &namebucket);

string reversecompletment(const string& str );

bool adjacent (const string& node1,const  string& node2,int k);

string readn(ifstream *file,uint64_t n);

int chartoint(char c);

string minimalsub(const string &w, const int &p,const int &k);

string minimalsub2(const string &w, const int &p,const int &k);

class neighbour
{
	public:
		array<pair<uint64_t,unsigned char>,8> list;
		//~ neighbour()
		//~ {
			//~ for(int )
			//~ list[i]=make_pair(0,0);
		//~ }
		uint64_t nbtype(unsigned char c);
		uint64_t gtype(unsigned char c);

		void add(uint64_t p,unsigned char b);
		unsigned char remove(uint64_t v);
		unsigned char removep(uint64_t v,unsigned char c);
		unsigned char removetype(unsigned char c);

};

class graph
{
	public:
		uint64_t n;
		int k;
		vector<string> nodes;
		vector<int> leftmins;
		vector<int> rightmins;
		unordered_multimap<uint64_t,uint64_t> map;
		unordered_multimap<uint64_t,uint64_t> maprev;
		vector<neighbour> neighbor;

		graph(const int ni)
		{
			k=ni;
			n=1;
			nodes.push_back("");
			leftmins.push_back(-1);
			rightmins.push_back(-1);
		}

		uint64_t getkey(string str);
		uint64_t getkeyrevc(string str);
		uint64_t becompacted(uint64_t nodeindice, int min, unsigned char *);
		int weight();
		void addvertex(const string str);
        void addleftmin(int mini);
        void addrightmin(int mini);
		void debruijn();
		void compressh(int min=-1);
		void compress();
		void importg(const char *name);
		void print(const char *name);
		void printedges(const char *name);
		void compact(uint64_t nodeindice,uint64_t with, unsigned char type);
		void reverse(int64_t with);
		void look(const uint64_t nodeindice, const string min);

};

#endif
