#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include "lm.h"
#include "ograph.h"
#include "debug.h"

using namespace std;

// this code crashes now; this is the old behavior when no arguments was specified, it runs some test
void test()
{
    srand(time(NULL));
    bool b(true);
    int n(1000);
    while(b)
    {

        int sys(0);
        int k(20);
        int m(2);

        cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
        createinputlm(n,k,"randomgenome");
        sys+=system("cat randomgenome | sort | uniq > input.dot");
        cout<<"GO!"<<endl;
        graph G(k),G2(k);
        G.importg("input.dot");
        G.debruijn();
        G.compress();
        G.print("output.dot");
        cout<<"GO!"<<endl;
        createoutfile("input.dot","outputlowmemory.dot",k,m);
        if(!checkfile("output.dot","outputlowmemory.dot",k)){
            cout<<"Errors occurred !"<<endl;
            b=false;
        }
        else{
            cout<<"Success !"<<endl;
            n*=10;
            cout<<n<<endl;
        }
        b=false;
    }

}



int main(int argc, char ** argv)
{
	int sys(0);
	if(argc==1)
	{
        printf("usage: <input> [output.dot] [minimizer length]\n");
        printf("Note: default behavior (minimizer length = 10) requires that you type 'ulimit -n 1100' in your shell prior to running bcalm, else the software will crash\n");
        exit(1);
	}
	if(argc==2)
	{
		string input(argv[1]);
		string output("compacted.dot");
		int m(5);
		if(testulimit(1100))
		{
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}
	if(argc==3)
	{
		string input(argv[1]);
		string output(argv[2]);
		int m(5);
		if(testulimit(1100))
		{
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}
	if(argc==4)
	{
		string input(argv[1]);
		string output(argv[2]);
		int m(atoi(argv[3])/2);
		if(testulimit(pow(4,m)+50))
		{
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}


	return sys;
}
