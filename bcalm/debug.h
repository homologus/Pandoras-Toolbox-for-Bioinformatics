#ifndef DEBUG
#define DEBUG

using namespace std;

void fastatodot(const char * name,const char * namout);

void createinputlm(int64_t lr,int k,const char *name);

bool checkfile(string name1, string name2,int k);

int detectk(const string& input);

bool testulimit(int l);

#endif
