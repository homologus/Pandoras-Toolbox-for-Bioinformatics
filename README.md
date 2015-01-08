<img src="https://cloud.githubusercontent.com/assets/3611756/5643160/a8574d64-9609-11e4-9abb-f5ee60779fa0.png" width="400" align="center">

Pandora's Toolbox
=================

License – GNU GPL3

Condition for use - Please do NOT cite Pandora's Toolbox. Instead cite papers from authors, whose programs you are using. Those links are given below.


Description
=================

Pandora’s Toolbox is a collection of source codes of well-known bioinformatics programs all at the same place. You can now run ‘git clone’ and then ‘make’ once in the main folder. That should compile everything (not there yet) and collect executables in the ‘bin’ folder.

Previously, every time we tried to build them, we faced the problem of having to collect the source codes from various websites. Then each program had to be compiled individually, and some of those using boost were nightmares to compile. Also, many programs have interdependencies and we ended up having five copies of different versions of BWA or samtools codes.

Our collection reduces those dependencies and makes life easier for those working on different bioinformatics programs. The collection includes a boost_1_55_0 folder to remove external boost dependencies. We are working on sorting out all bwa and other inter-dependencies.

The following programs are included in the current version. Please DO NOT cite us, but cite the authors of individual programs.

1. KMC2
========

KMC2 is an efficient kmer counter that uses minimizers.

Cite -

<A href=http://arxiv.org/abs/1407.1507>KMC 2: Fast and resource-frugal k-mer counting - Sebastian Deorowicz, Marek Kokot, Szymon Grabowski, Agnieszka Debudaj-Grabysz</A>

http://sun.aei.polsl.pl/kmc/
http://sun.aei.polsl.pl/kmc/download/kmc.tar.gz

**** License *****
* KMC software distributed under GNU GPL 2 licence.

* libbzip2 is open-source (BSD-style license)

* gzip is free, open-source

* asmlib is under the licence GNU GPL 3 or higher

* Boost is under Boost Software License (free, open-source), see http://www.boost.org/users/license.html.

In case of doubt, please consult the original documentations.


2. bwa (included BWAMEM)
============================

BWA is a short read aligner.

Cite -
<A href=http://arxiv.org/abs/1303.3997>Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM - Heng Li</A>

https://github.com/lh3/bwa


3. DALIGNER
============

DALIGNER is an aligner for long noisy reads.

Cite Gene Myers, WABI 2014.

https://github.com/thegenemyers

DALIGNER - Find all significant local alignments between reads
DEXTRACTOR -
DAZZ_DB	- The Dazzler Data Base


4. HMMER3
==========

Cite  –

<A href=http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002195>Accelerated Profile HMM Searches - Sean R. Eddy</A>



5. bcalm – de Bruijn graph compressor
=====================================

This code compresses kmers into dBG using Minimizers. BCALM's license is BSD with Attribution. 
In informal words, any use is permitted as long as it acknowledges the original authors.

Cite -

<A href=http://arxiv.org/abs/1401.5383>On the representation of de Bruijn graphs - Rayan Chikhi, Antoine Limasset, Shaun Jackman, Jared Simpson, Paul Medvedev</A>

https://github.com/Malfoy/bcalm



6. Minia – assembler
=====================

inia is a short-read assembler based on a de Bruijn graph, capable of assembling a human genome on a desktop computer in a day. The output of Minia is a set of contigs. Minia produces results of similar contiguity and accuracy to other de Bruijn assemblers (e.g. Velvet).


Cite -

R. Chikhi, G. Rizk. Space-efficient and exact de Bruijn graph representation based on a Bloom filter, WABI 2012

@inproceedings{minia,
 author = {Chikhi, Rayan and Rizk, Guillaume},
 title = {Space-Efficient and Exact de Bruijn Graph Representation Based on a Bloom Filter.},
 booktitle = {WABI},
 pages = {236-248},
 publisher = {Springer},
 series = {Lecture Notes in Computer Science},
 volume = 7534,
 year = 2012
} 

https://gatb.inria.fr/



7. samtools – NGS data handling
===============================

Cite Heng Li.

https://github.com/samtools


8. SOAPdenovo2 – genome assembler
=================================

Cite Luo et al.

GPL3.0 license.

http://softlayer-dal.dl.sourceforge.net/project/soapdenovotrans/SOAPdenovo-Trans/src/v1.04/SOAPdenovo-Trans-src-v1.04.zip

9. SOAPdenovo-trans – transcriptome assembler
=============================================

Cite Xie et al.


10. SPAdes – genome assembler
=============================

Cite the SPAdes team.

http://bioinf.spbau.ru/en/spades

http://bioinf.spbau.ru/en/content/spades-download-0


11. Trinity – transcriptome assembler
=====================================

Cite Grabher et al.

http://trinityrnaseq.sourceforge.net/

http://sourceforge.net/projects/trinityrnaseq/files/
http://hivelocity.dl.sourceforge.net/project/trinityrnaseq/trinityrnaseq_r20140717.tar.gz


12. sailfish – RNAseq expression analysis
==========================================

Cite Rob Patro et al..
[We are still working on compilation of this code.]

http://www.cs.cmu.edu/~ckingsf/software/sailfish/

https://github.com/kingsfordgroup/sailfish


13. RAPSearch2
===============


14. Tophat
============


15. Cufflinks
==============




About the name
==============

We named the modules after Symbion pandora, an unusual animal discovered by Danish researchers Peter Funch and R. M. Christensen in 1995. Many thanks to professor Funch for sharing the SEM image at the top of this post.

More on Symbion pandora –

"Symbion pandora is a jug-shaped microscopic aquatic animal that dwells on the mouth-parts of Norway lobsters. The animals are less than ½ mm wide, with sac-like bodies, and three distinctly different forms in different parts of their three-stage life cycle.

They are so unlike any known animal that its discovery by Danish scientists in 1995[1] led to the creation of a new phylum. The phylum Cycliophora, from the Greek for ‘carrying a small wheel’, was named after the creature’s circular mouth."

From their original Nature paper –

"THE mouthparts of the Norway lobster Nephrops are colonized by an acoelomate metazoan, Symbion pandora gen. et sp. nov. Sessile stages continually produce inner buds replacing feeding structures. They also produce one of three motile stages: (1) larvae containing new feeding stages, (2) dwarf males, which settle on feeding stages, or (3) females, which settle onto lobster mouthparts, and eventually degenerate, giving rise to dispersive larvae. All motile stages are short-lived, and do not feed. The structure and function of the cilia suggest a phylogenetic position in Protostomia, while some aspects of inner budding and brooding of larvae are similar to those of Entoprocta and Ectoprocta. The dispersive larva possesses a mesodermal supporting chordoid structure, otherwise absent in protostomian larvae. We believe that all the above features of this previously undescribed species warrant the recognition of a new phylum with affinities to Ectoprocta and Entoprocta."
