#ifndef INPUT
#define INPUT

#include <string>



using namespace std;

void init_input(string namefile, const int k);

string next_input(const int k);


class InputDot{
	
uint64_t size, buffer_size, nb_buffers, index,indexb;
string buffer,kmer;
ifstream in;

public:
void init_input(string namefile, const int k);

string next_input(const int k);

};

#endif
