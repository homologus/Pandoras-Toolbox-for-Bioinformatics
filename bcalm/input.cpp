#include <string>
#include <iostream>
#include <fstream>

#include "input.h"
#include "ograph.h"



using namespace std;




void InputDot::init_input(string namefile, const int k){
	in.open(namefile);
	in.seekg(0,ios_base::end);
	size = in.tellg();
    buffer_size = 10*(k+2);
    nb_buffers = size/buffer_size;
	in.seekg(0,ios::beg);
    index = 0;
    indexb=0;
    if(nb_buffers!=0)
    {
		buffer=readn(&in,buffer_size);
	}
	else
	{
		buffer=readn(&in,size-nb_buffers*buffer_size);
	}
}

string InputDot::next_input(const int k){
	kmer="";
	if(index>=buffer.size())
	{
		index=0;
		indexb++;
		if(indexb<nb_buffers)
		{
			buffer=readn(&in,buffer_size);
		}
		else
		{
			if(indexb==nb_buffers)
			{
				buffer=readn(&in,size-nb_buffers*buffer_size);
			}
			else{
				buffer="";
			}
		}
	}

	if(buffer.size()>=index+k)
	{
		kmer=buffer.substr(index,k);
		index += k+2;
	}
	return kmer;
	
}

// l'idee sera d'avoir des classes ayant la meme intreface que InputDot pour gerer d'autres format d'entree
