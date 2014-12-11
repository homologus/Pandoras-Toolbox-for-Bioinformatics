#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>

/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

using namespace std;


static inline int nt2num(char c){
    // 'a' -> 0, 'c' -> 1; 'g' -> 2, 't' -> 3
    // inspired by G. Rizk
    char d = (c >> 1) & 3;
    if (d > 1) 
        d ^= 1;
    return d;
}

char num2nt(int num){
	if (num == 0)
			return 'a';
	else if (num == 1)
		return 'c';
	else if (num == 2)
		return 'g';
	else if (num == 3)
		return 't';
	assert(0);
	return '*';
}

static inline int nt2num(char c, int pos){
	int offset = pos % 4;
	return (nt2num(c) + offset ) % 4;
}


static inline char num2nt(int num, int pos){
	int offset = pos % 4;
	return num2nt((num - offset + 4) % 4);
}

int shash(const string& s, int& previous_shash, unsigned int start_pos, int length){
    if (length == -1)
    {
        length = s.length();
    }
	int retval = 0, offset = 0;
    if (previous_shash != -1)
    {
        // shortcut, assuming consecutive m-mers
        retval = (previous_shash >> 2);
        offset = length - 1;
    }
	for (unsigned int pos = start_pos; pos+ offset < start_pos + length; pos++){
		//int num = nt2num(s[pos], pos - start_pos); // would need to be adapted to previous_hash
		int num = nt2num(s[pos + offset], 0);
		//retval = retval + pow(4, s.length() - pos - 1) * num;
        retval +=  (num  << (2*(pos + offset - start_pos))); // actually compute the reverse 2bit representation now
	}
    previous_shash = retval;

    // debug
    // printf("shash %s -> %0X\n",s.substr(start_pos,length).c_str(),retval);
	return retval;
}

string inverse_shash(int num, int len){
	string s(len, 'X');
    int original_num = num;
	for (int pos = 0; pos < len; pos++){
		int val = num % 4;
		num = (num - val) / 4;
		s[pos] = num2nt(val);
	}

    // debug
    // printf("inverse shash %0X -> %s\n", original_num, s.c_str());
	return s;
}

uint32_t *m_mer_counts;

void init_m_mers_table(int m)
{
    m_mer_counts = new uint32_t[(int)pow(4,m)];
    //printf("Allocating %d MB for m-mers\n",(((int)pow(4,m)*sizeof(uint32_t))/1024)/1024);
    for (int i = 0; i < (int)pow(4,m); i++)
        m_mer_counts[i] = 0;
}

void count_m_mers(string str, int m, int k)
{
    int previous_shash = -1;
	for(int i(0) ;i < k-m+1 ; i+=1){
        m_mer_counts[shash(str, previous_shash, i, m)] ++;
    }
}

uint32_t *best_possible_hash_function;

void create_hash_function_from_m_mers(int m)
{
    int rg = pow(4,m);
    vector<pair<int, int> > counts;
    for (int i(0); i < rg; i++)
    {
        if (m_mer_counts[i] > 0)
            counts.push_back(make_pair(m_mer_counts[i],i));
    }

    printf("Sorting m-mer counts\n");
    sort(counts.begin(),counts.end());

    cout<< counts[counts.size()-1].first << " " << inverse_shash(counts[counts.size()-1].second,m) <<endl;
    cout<< counts[counts.size()-2].first << " " << inverse_shash(counts[counts.size()-2].second,m) <<endl;
    cout<< counts[counts.size()-3].first << " " << inverse_shash(counts[counts.size()-3].second,m) <<endl;

    delete[] m_mer_counts; // fixme: restore this
    best_possible_hash_function = new uint32_t[rg];

    for (int i = 0; i < rg ; i++)
        best_possible_hash_function[i] = 0;

    for (unsigned int i = 0; i < counts.size(); i++)
        best_possible_hash_function[counts[i].second] =  i;
}

// not needed anymore
HashMap build_hash_map(int len){
	HashMap hm;
    return hm; // deactivate it

	for(int i(0);i<pow(4,len); i++){
		string s=inverse_shash(i,len);
		hm[s]=i;
	}
	hm.rehash(2 *(pow(4,len)));
	return hm;
}



int getash(const string& s, int& previous_shash, int start_pos = 0, int length = -1){
	//return ((*hm)[s]); // slow
    //return shash(s, start_pos, length);
    return best_possible_hash_function[shash(s, previous_shash, start_pos, length)];
}

int minimiserv(const string &node,const int &minimisersize){
    int previous_shash = -1;
	int minimiser_value(getash(node, previous_shash, 0, minimisersize)),vsub;
	for(uint64_t i(1);i<node.size()-minimisersize+1;i++){
		vsub=getash(node, previous_shash,i,minimisersize);
		if( minimiser_value > vsub){
			minimiser_value = vsub;
		}
	}
    return(minimiser_value);
}

int minimiserrc(const string &node,const int &minimisersize){
	int h1, h2;
	h1 = minimiserv(node,minimisersize);
	h2 = minimiserv(reversecompletment(node),minimisersize);
    return (h1 > h2) ? h2 : h1;
}

int minimiserrc_openmp(const string &node,const int &minimisersize){
	string nodes[2];
	int resHash[2];

	nodes[0] = node;
	nodes[1] = reversecompletment(node);

	//~ #pragma omp parallel for
	for (int i=0; i<2; i++){
		resHash[i] = minimiserv(nodes[i],minimisersize);
	}

	if(resHash[0] < resHash[1]){
		return resHash[0];
	}
	else{
		return resHash[1];
	}
}

int minbutbiggerthan(int m1, int m2, const string &namebucket){
	int h(stoi(namebucket));
	if(m1<m2){
		if(m1>h)
			return m1;
		if(m2>h)
			return m2;
	}
	else{
		if(m2>h)
			return m2;
		if(m1>h)
			return m1;
	}
	return -1;
}


string reversecompletment(const string& str)
{
	string res(str);
    int n = str.size();
	for(int i(n-1), j(0); i > -1; i--, j++)
    {
        unsigned char c = str[i];
        unsigned char d = (c >> 4)&7;
        if (d >= 6) //(c >= 'a' && c <= 't')
        {
            // translates acgt to tgca
            c ^= 4;
            if ((c&3) != 3)
                c ^= 17;
            res[j] = c;
            continue;
        }
        if (d == 2)  // switch '+' with '-'
            res[j] = c ^ 6;
        else
        {
            // else it will only be a number, just copy it
            res[j] = c;
        }
    }
    return res;

}

//read n character
string readn(ifstream *file,uint64_t n){
	string contents(n,'0');
	file->read(&contents[0],n);
	return(contents);
}



bool adjacent(const string& node1,const  string& node2,int k){
	return(node1.substr(node1.size()-k+1,k-1)==node2.substr(0,k-1));
}



int chartoint(char c){
	switch(c){
		case 'a':
		return 0;
		case 'c':
		return 1;
		case 'g':
		return 2;
		case 't':
		return 3;
		default:
		cout<<"Problem with chartoint:"<<c<<endl;
		assert(0);
		return 0;
	}
}

uint64_t stringtoint(const string& str){
	uint64_t res(0);
	for(uint64_t i(0);i<str.size();i++){
		res<<=2;
		res+=chartoint(str[i]);
	}
	return res;
}

uint64_t stringtointc(const string& str){
	uint64_t res(0);
	for(int64_t i(str.size()-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


void graph::addvertex(string str){
	n++;
	nodes.push_back(str);
	uint64_t i(nodes.size()-1);
	uint64_t key(getkey(str));
	uint64_t keyrc(getkeyrevc(str));
	map.insert({key,i});
	maprev.insert({keyrc,i});
}

void graph::addleftmin(int mini){
	leftmins.push_back(mini);
}

void graph::addrightmin(int mini){
	rightmins.push_back(mini);
}

void graph::debruijn(){
	neighbor=vector<neighbour> (n);
	string node,kmmer,kmmerr;
	uint64_t key,keyrc;
	for(uint64_t i(1);i<n;i++){
		node=nodes[i];
		key=stringtoint(node.substr(node.size()-k+1,k-1));
		keyrc=stringtointc(node.substr(0,k-1));
		auto it(map.equal_range(key));
		for(auto j(it.first);j!=it.second;j++){
			//if k>32 collision can occur
			if(adjacent(node,nodes[j->second],k)){
				neighbor[i].add(j->second,1);
				neighbor[j->second].add(i,4);
			}
		}
		it=(maprev.equal_range(key));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(node,reversecompletment(nodes[j->second]),k)){
				neighbor[i].add(j->second,2);
				neighbor[j->second].add(i,2);
			}
		it=(map.equal_range(keyrc));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(reversecompletment(node),nodes[j->second],k)){
				neighbor[i].add(j->second,3);
				neighbor[j->second].add(i,3);
			}
	}
	map.clear();
}


bool accordtomin(int min, int left_or_right_min){
	if(min == -1){
		return true;
	}

    if(left_or_right_min==min)
		return true;

	return false;

}

uint64_t graph::becompacted(uint64_t nodeindice, int min, unsigned char *type){

	*type=0;
	string node=nodes[nodeindice];

	if(node.empty())
		return 0;

    int leftmin = leftmins[nodeindice];
	int rightmin = rightmins[nodeindice];

	auto neigh(neighbor[nodeindice]);
	int one(neigh.nbtype(1)),two(neigh.nbtype(2)),three(neigh.nbtype(3)),four(neigh.nbtype(4));
	int in(three+four),out(one+two);
	if(out==1 && accordtomin(min,rightmin)){
		if(one==1){
			uint64_t sonindice(neigh.gtype(1));
			*type=1;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		else{
			uint64_t sonindice(neigh.gtype(2));
			*type=2;
			if(neighbor[sonindice].nbtype(2)+neighbor[sonindice].nbtype(1)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	if(in==1 && accordtomin(min,leftmin)){
		if(three==1){
			uint64_t sonindice(neigh.gtype(3));
			*type=3;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		else{
			uint64_t sonindice(neigh.gtype(4));
			*type=4;
			if(neighbor[sonindice].nbtype(1)+neighbor[sonindice].nbtype(2)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	return 0;
}

void graph::reverse(int64_t with){
	string newnode(nodes[with]);
	int64_t indice,type;
	if(newnode>(reversecompletment(newnode))){
		nodes[with]=reversecompletment(newnode);

        int temp = leftmins[with];
        leftmins[with] = rightmins[with];
        rightmins[with] = temp;

		for(auto i(0);i<8;i++){
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;

			if(type==1){
				neighbor[indice].removep(with,4);
				neighbor[with].removep(indice,1);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
				continue;
			}
			if(type==2){
				neighbor[indice].removep(with,2);
				neighbor[with].removep(indice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
				continue;
			}
			if(type==3){
				neighbor[indice].removep(with,3);
				neighbor[with].removep(indice,3);
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
				continue;
			}
			if(type==4){
				neighbor[indice].removep(with,1);
				neighbor[with].removep(indice,4);
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
				continue;
			}
		}
	}
}


void graph::compact(uint64_t nodeindice,uint64_t with, unsigned char c){
	string newnode,node(nodes[nodeindice]),son(nodes[with]);
	uint64_t indice;
	if(nodeindice==with || node.empty() || son.empty())
		return;
	unsigned char type;
	switch(c){
		case 1:
		newnode=node+son.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;

        leftmins[with] = leftmins[nodeindice];
        //rightmins[with] = rightmins[with];

		for(int  i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(indice==nodeindice){
				continue;
			}
			if(type==3 ){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 2:
		newnode=node+reversecompletment(son).substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;

        rightmins[with] = leftmins[with];
        leftmins[with] = leftmins[nodeindice];

		for(auto i(0);i<8;i++){
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;
			neighbor[indice].remove(with);
			neighbor[with].remove(indice);
			if(type==3){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==4){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(type==3){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 3:
		newnode=reversecompletment(node)+son.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;

        leftmins[with] = rightmins[nodeindice];
        //rightmins[with] = rightmins[with];


		neighbor[with].removep(nodeindice,3);
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			if(type==1){
				neighbor[indice].removep(nodeindice,4);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==2 ){
				neighbor[indice].removep(nodeindice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 4:
		newnode=son+node.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;

        // leftmins[with] = leftmins[with];
        rightmins[with] = rightmins[nodeindice];

		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(type==1){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==2 ){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
		break;
	}
}


//Compact the graph but not the nodes that should be compacted in an other bucket
void graph::compressh(int min){
	unsigned char type(0);
	uint64_t with;
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++){
		with=becompacted(nodeindice,min,&type);
		if(with!=0)
			compact(nodeindice,with,type);
	}
}

//Compact the graph but not the nodes that should be compacted in an other bucket
void graph::compress(){
	HashMap hm;
	unsigned char type(0);
	uint64_t with;
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++){
		with=becompacted(nodeindice,-1,&type);
		if(with!=0)
			compact(nodeindice,with,type);
	}
}



//import graph from file
void graph::importg(const char *name){
	cout<<"importg"<<endl;
	ifstream fichier(name,ios::in);
	string line,vertex;
	if(fichier)
		while(!fichier.eof()){
			getline(fichier,line);
			int64_t lp(0);
			for(unsigned int i(0);i<line.size();i++){
				if(line[i]==';'){
					addvertex(line.substr(lp,i-lp));
					lp=i+1;
				}
			}
		}
	else
		cerr<<"no file"<<endl;
}



int graph::weight(){
	int res(0);
	for(auto k(nodes.begin());k!=nodes.end();k++)
		res+=k->size();
	return res;
}



void graph::print(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier)
		for(uint64_t i(1);i<n;i++){
			string s(nodes[i]);
			if(s!="")
				fichier<<s<<";"<<endl;
		}
}

void graph::printedges(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier){
		fichier << "digraph test {" <<endl;
		for(uint64_t i(1);i<n;i++){
			string s(nodes[i]);
			if(s!=""){
				fichier<<s<<";"<<endl;
				//~ auto v=neighbor[i].son;
				//~ for(auto j=v.begin();j!=v.end();j++)
					//~ if(j->first!=0)
						//~ if(nodes[j->first]!="")
							//~ fichier<<s<<"->"<<nodes[j->first]<<";"<<endl;
			}
		}
		fichier << "}"<<endl;
		fichier.close();
	}
}


uint64_t graph::getkey(string str){
	return stringtoint(str.substr(0,k-1));
}


uint64_t graph::getkeyrevc(string str){
	return stringtointc(str.substr(str.size()-k+1,k-1));
}


void neighbour::add(uint64_t p,unsigned char c){
	for(int i(0);i<8;i++){
		if(list[i].first==0  ){
			list[i]=make_pair(p,c);
			return;
		}
		if(list[i].first==p && list[i].second==c )
			return;
	}
}



uint64_t neighbour::nbtype(unsigned char c){
	uint64_t ret(0);
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			ret++;
	return ret;
}




uint64_t neighbour::gtype(unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			return list[i].first;
	cout<<"Bug with neighbour"<<endl;
	return 0;
}




unsigned char neighbour::remove(uint64_t v){
	for(int i(0);i<8;i++)
		if(list[i].first==v){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	return 0;
}

unsigned char neighbour::removep(uint64_t v,unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].first==v && list[i].second==c){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	return 0;
}


unsigned char neighbour::removetype(unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			list[i].first=0;
	return 0;
}

