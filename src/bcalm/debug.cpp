#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "debug.h"
#include "ograph.h"



using namespace std;


void fastatodot(const char * name,const char * namout){
	ifstream in(name);
	ofstream out(namout);
	string data;
	while(!in.eof()){
		getline(in,data);
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		if(data!="")
			out<<data<<";"<<endl;
	}
}

int detectk(const string& input){
	ifstream in(input);
	string line;
	if(in){
		getline(in,line);
	}
	return(line.size()-1);
}

bool testulimit(int l)
{
	system("ulimit -n > .ulimit");
	ifstream in(".ulimit");
	string line;
	getline(in,line);
	int n;
	n=stoi(line);
	remove(".ulimit");
	return (n>=l);
}


void createinputlm(int64_t lr,int k,const char *name){
	ofstream out(name,ios::trunc);
	int r;
	string c;
	string kmer(k,'a');
	for(int b(0);b<k;b++){
		r=rand()%4;
		switch(r){
			case 1:{
				kmer[b]='a';
				break;
			}
			case 2:{
				kmer[b]='c';
				break;
			}
			case 3:{
				kmer[b]='g';
				break;
			}
			case 0:{
				kmer[b]='t';
				break;
			}
		}
	}
	for(int64_t b(0);b<lr;b++){
		kmer=kmer.substr(1,k-1);
		r=rand()%4;
		switch(r){
			case 1:{
				c='a';
				break;
			}
			case 2:{
				c='c';
				break;
			}
			case 3:{
				c='g';
				break;
			}
			case 0:{
				c='t';
				break;
			}
		}
		kmer+=c;
		if(kmer<reversecompletment(kmer)){
			out<<kmer<<";"<<endl;
		}
		else{
			if(kmer>reversecompletment(kmer)){
				out<<reversecompletment(kmer)<<";"<<endl;
			}
		}
	}
}


bool checkfile(string name1, string name2,int k){
	int fail(0),lowfail(0);
	ifstream t1(name1), t2(name2);
	string line,kmer;
	unordered_map<string,bool> s1,s2,e1,e2;
	while(!t1.eof()){
		getline(t1,line);
		if(line.size()>2){
			string node=line.substr(0,line.size()-1);
			node=min(node,reversecompletment(node));
			s1.insert(make_pair(node,false));
		}
	}
	while(!t2.eof()){
		getline(t2,line);
		if(line.size()>2){
			string node=line.substr(0,line.size()-1);
			node=min(node,reversecompletment(node));
			s2.insert(make_pair(node,false));
		}
	}
	for(auto it=s1.begin(); it!=s1.end(); it++){
			string str=it->first;
			auto smt=s2.find(str);
			if(smt==s2.end()){
				lowfail++;
				for(int i(0);i+k<=(int)str.size();i++){
					kmer=str.substr(i,k);
					kmer=min(kmer,reversecompletment(kmer));
					e2.insert(make_pair(kmer,false));
				}
			}
			else
				smt->second=true;
	}
	for(auto it=s2.begin(); it!=s2.end(); it++)
		if(!it->second)
		{
			string str=it->first;
			lowfail++;
			for(int i(0);i+k<=(int)str.size();i++){
				kmer=str.substr(i,k);
				kmer=min(kmer,reversecompletment(kmer));
				e1.insert(make_pair(kmer,false));
			}
		}
	for(auto it=e2.begin(); it!=e2.end(); it++){
		string str=it->first;
		auto smt=e1.find(str);
		if(smt==s1.end()){
			fail++;
			cout<<str<<endl;
		}
		else
			smt->second=true;
	}
	for(auto it=e1.begin(); it!=e1.end(); it++){
		if(!it->second){
			fail++;
			cout<<it->first<<endl;
		}
	}
	cout<<"Lowrrors:"<<lowfail<<endl;
	cout<<"Errors:"<<fail<<endl;
	return(lowfail==0);
}
