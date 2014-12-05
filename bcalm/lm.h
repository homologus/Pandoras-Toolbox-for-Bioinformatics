#ifndef LM
#define LM
#include "ograph.h"

using namespace std;

void createoutfile(const char *namein,const char *nameout,const int k,const int m);

void sortentry(string namefile, const int k, const int m, bool create_buckets = true, bool m_mer_count = false);

#endif
