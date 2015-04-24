<img src="https://cloud.githubusercontent.com/assets/3611756/5643160/a8574d64-9609-11e4-9abb-f5ee60779fa0.png" width="400" align="center">

Pandora's Toolbox
=================

License – GNU GPL3


Condition for use
=================

Please do NOT cite Pandora's Toolbox. Instead cite papers from individual authors, whose programs you are using. Those links are given below.


Description
=================

Pandora’s Toolbox is a collection of source codes of well-known bioinformatics programs all at the same place. You can now run ‘git clone’ and then ‘make’ once in the main folder. That should compile everything (not there yet) and collect executables in the ‘bin’ folder.

Previously, every time we tried to build them, we faced the problem of having to collect the source codes from various websites. Then each program had to be compiled individually, and some of those using boost were nightmares to compile. Also, many programs have interdependencies and we ended up having five copies of different versions of BWA or samtools codes.

Our collection reduces those dependencies and makes life easier for those working on different bioinformatics programs. The collection includes a boost_1_55_0 folder to remove external boost dependencies. We are working on sorting out all bwa and other inter-dependencies.

The following programs are included in the current version. Please DO NOT cite us, but cite the authors of individual programs.

1. KMC2
========

KMC2 is an efficient kmer-counter that does not require significant RAM. It is disk-based and uses a minimizer-like method to partition the reads.

1.  Deorowicz, S., Kokot, M., Grabowski, Sz., Debudaj-Grabysz, A., KMC 2: Fast and resource-frugal k-mer counting, Bioinformatics, 2015.

2. Deorowicz, S., Debudaj-Grabysz, A., Grabowski, Sz., Disk-based k-mer counting on a PC, BMC Bioinformatics, 2013; 14():Article no. 160.


2. BWA and BWA-MEM
============================

BWA-MEM searches for given short and long reads within an existing sequence or collection of sequences. For example, it can be used to find matches of millions of short reads in the human genome.

1. [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. Heng Li (2013).](http://arxiv.org/abs/1303.3997)


3. DALIGNER
============

DALIGNER finds overlaps between long noisy reads with extensive amount of insertion-deletion errors.

Eugene Myers, WABI 2014 (Workshop on Algorithms in Bioinformatics), Sept 8-10, Wroclaw, Poland.


4. HMMER3
==========

HMMER is a sensitive nucleotide and protein search program. It uses hidden Markov model.

1. Accelerated profile HMM searches.  S. R. Eddy.  PLoS Comp. Biol., 7:e1002195, 2011.


5. bcalm – de Bruijn graph compressor
=====================================

Splitting the short reads into k-mers and building the de Bruijn graph structure is usually the first step of assembly for de Bruijn graph-based algorithms. This step is often very memory intensive. A compressed de Bruijn graph combines contiguous k-mers into longer sequences and can be held with less memory. BCALM is a memory-efficient program to generate compressed de Bruijn graph from a collection of k-mers.

1.	[On the representation of de Bruijn graphs, Rayan Chikhi et al. (2014)](http://arxiv.org/abs/1401.5383).

6. Minia – assembler
=====================

Minia is a contig assembler for short reads. It uses very low amount of memory, such as less than 4GB to assemble the human genome.

1. R. Chikhi, G. Rizk. Space-efficient and exact de Bruijn graph representation based on a Bloom filter, WABI 2012.

2. K. Salikhov, G. Sacomoto and G. Kucherov. Using cascading Bloom filters to improve the memory usage for de Brujin graphs, WABI 2013.

7. samtools – NGS data handling
===============================

Samtools is useful for processing short-read alignment data.

https://github.com/samtools

8. SOAPdenovo2 – genome assembler
=================================

SOAPdenovo uses de Bruijn graph-based algorithms to assemble large eukaryotic genomes from short read libraries.

1. [De novo assembly of human genomes with massively parallel short read sequencing. R. Li et al. Genome Res 20, 265-72 (2010).](http://www.ncbi.nlm.nih.gov/pubmed/20019144)

2. [SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler R. Luo et al. GIGAscience (2012).](http://www.gigasciencejournal.com/content/1/1/18)


9. SOAPdenovo-trans – transcriptome assembler
=============================================

SOAPdenovo-trans assembles short reads of RNAseq experiments into transcripts.

1. [SOAPdenovo-Trans: De novo transcriptome assembly with short RNA-Seq reads, Yinlong Xie et al. (2014).](http://arxiv.org/ftp/arxiv/papers/1305/1305.6760.pdf)

10. SPAdes – genome assembler
=============================

SPAdes, developed by Algorithmic Biology lab in Saint Petersburg, Russsia is among the most efficient and versatile genome assemblers that uses de Bruijn graphs. SPAdes was originally designed to assemble single-celled bacterial genomes, but it appears to work well for multi-cell data, as well as small to mid-sized eukaryotic genomes.

1. Anton Bankevich, Sergey Nurk, Dmitry Antipov, Alexey A. Gurevich, Mikhail Dvorkin, Alexander S. Kulikov, Valery M. Lesin, Sergey I. Nikolenko, Son Pham, Andrey D. Prjibelski, Alexey V. Pyshkin, Alexander V. Sirotkin, Nikolay Vyahhi, Glenn Tesler, Max A. Alekseyev, and Pavel A. Pevzner. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19(5) (2012), 455-477. doi:10.1089/cmb.2012.0021


11. Trinity – transcriptome assembler
=====================================

Trinity is a stand alone *de novo* transcriptome assembler that uses de Bruijn graph-based algorithms.

1.  [Full-length transcriptome assembly from RNA-Seq data without a reference genome, Manfred G Grabherr et al. Nature Biotechnology 29, 644.652 (2011)](http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html)


12. sailfish – RNAseq expression analysis
==========================================

Sailfish is a lightweight program for quantifying the abundance of previously annotated RNA isoforms in RNAseq data.

1. [Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms,  Rob Patro, Stephen M Mount	& Carl Kingsford, Nature Biotechnology 32, 462-464 (2014).](http://www.nature.com/nbt/journal/v32/n5/abs/nbt.2862.html) 


13. RAPSearch2
===============

RAPsearch2 searches for matches to a protein sequences in a database of sequences. It is 80x faster than BLAST without significant loss of quality.

1. [Yuzhen Ye, Jeong-Hyeon Choi and Haixu Tang. RAPSearch: a Fast Protein Similarity Search Tool for Short Reads. BMC Bioinformatics 2011, 12:159.](http://www.biomedcentral.com/1471-2105/12/159/abstract).

2. [RAPSearch2: a fast and memory-efficient protein similarity search tool for next-generation sequencing data Yongan Zhao, Haixu Tang and Yuzhen Ye. 2011.](http://bioinformatics.oxfordjournals.org/content/28/1/125.full)


14. Tophat and Cufflinks
============

Tophat and Cufflinks process RNAseq reads aligned onto a reference genome and resolve the intron-exon junctions.

1.	[TopHat: discovering splice junctions with RNA-Seq Cole Trapnell, Lior Pachter and Steven L. Salzberg Bioinformatics. 2009](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672628/)

2.	Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks Cole Trapnell et al. Nature Protocols 7, 562.578 (2012).




About the name
==============

We named the modules after Symbion pandora, an unusual animal discovered by Danish researchers Peter Funch and R. M. Christensen in 1995. Many thanks to professor Funch for sharing the SEM image at the top of this post.

More on Symbion pandora –

"Symbion pandora is a jug-shaped microscopic aquatic animal that dwells on the mouth-parts of Norway lobsters. The animals are less than ½ mm wide, with sac-like bodies, and three distinctly different forms in different parts of their three-stage life cycle.

They are so unlike any known animal that its discovery by Danish scientists in 1995[1] led to the creation of a new phylum. The phylum Cycliophora, from the Greek for ‘carrying a small wheel’, was named after the creature’s circular mouth."

From their original Nature paper –

"THE mouthparts of the Norway lobster Nephrops are colonized by an acoelomate metazoan, Symbion pandora gen. et sp. nov. Sessile stages continually produce inner buds replacing feeding structures. They also produce one of three motile stages: (1) larvae containing new feeding stages, (2) dwarf males, which settle on feeding stages, or (3) females, which settle onto lobster mouthparts, and eventually degenerate, giving rise to dispersive larvae. All motile stages are short-lived, and do not feed. The structure and function of the cilia suggest a phylogenetic position in Protostomia, while some aspects of inner budding and brooding of larvae are similar to those of Entoprocta and Ectoprocta. The dispersive larva possesses a mesodermal supporting chordoid structure, otherwise absent in protostomian larvae. We believe that all the above features of this previously undescribed species warrant the recognition of a new phylum with affinities to Ectoprocta and Entoprocta."
