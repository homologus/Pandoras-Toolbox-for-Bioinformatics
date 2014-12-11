bcalm
=====

de Bruijn CompAction in Low Memory

Possible use :

    ./bcalm input.dot
compact input.dot with l=10 in compacted.dot

    ./bcalm input.dot output.dot
compact input.dot with l=10 in output.dot

    ./bcalm input.dot output.dot 8
compact input.dot with l=8 in output.dot



Nota Bene :   
Higher l mean lower memory but the algorithm will NOT work with l>10.   
Also note that to use l>=8 you have to allow 260 open files (type `ulimit -n 260` in a terminal)
and 1110 for l=10 (`ulimit -n 1100`).

Input
=====

Each line in the input file is a distinct k-mer in lower-case, ending with a semi-column. Here is a sample input.dot:

    actg;
    ctga;
    tgac;

You can convert DSK's output to bcalm's input as follows: use the `parse_results` program
provided with DSK with the `--dot` flag, e.g.: 

`./parse_results dsk_output.solid_kmers_binary --dot > input_bcalm.dot`

Output
=====

Each line is a simple path of the de Bruijn graph, in a similar format as the input.

License
=======

BCALM's license is BSD with Attribution. In informal words, any use is permitted as long as it acknowledges us. For example with a citation to:

    @inproceedings{reprdbg14,
    author = {Chikhi, R. and Limasset, A. and Jackman, S. and Simpson, J. and Medvedev, P.},
    title = {{On the representation of de Bruijn graphs}},
    booktitle = {RECOMB},
    publisher = {Springer},
    year = 2014,
    }

