# Complete-contigs
An assembler of unitigs, Y-to-V contigs and omnitigs.

# 1. Installation

In *src directory*, run

	make

This places the executable *complete-contigs* in the directory *bin*.

If you want to install the program for removing non-maximal omnitigs, In *src directory*, run

	make maximality

# 2. Usage

	complete-contigs OPTIONS

Brief example:

	complete-contigs -i <input .fa file>  -k 31 -g circular [-a 2] [-t 8]

Options:

	-h, --help          show this help message and exit
	--version           show program's version number and exit
	-i STRING, --input=STRING
						input .fa file with the single genome sequence
	-g STRING, --genome-type=STRING 
						genome type: linear|circular (default: circular)
	-k INT, --kmersize=INT
						kmer size (default: 31)
	-a INT, --abundance=INT
						minimum abundance for a k-mer to be present 
						in the de Bruijn graph (default: 1)
	-b, --build-only    build the de Bruijn graph and then exit
	-t INT, --threads=INT
    					number of threads, in [1..16]  (default: 1)

# 3. Removing non-maximal omnitigs

The outputted unitigs and Y-to-V contigs are maximal, but omnitigs are not necessarily. To remove non-maximal omnitigs, run for example
	
	bin/maximality chr10.k55.circular.omnitigs
