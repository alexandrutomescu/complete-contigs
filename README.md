# Complete-contigs
An assembler of unitigs, Y-to-V contigs and omnitigs

# Installation

In *src directory*, run

	make

This places the executable *complete-contigs* in the directory *bin*.

If you want to install the program for removing non-maximal omnitigs, In *src directory*, run

	make maximality

# Usage

	complete-contigs OPTIONS

Brief example:

	complete-contigs -i <input .fa file>  -k 31 -g circular [-a 2] [-t 8]

Options:

	-h, --help          show this help message and exit
	--version           show program's version number and exit
	-i STRING, --input=STRING
						input file, without extension (default assumed extension .fa)
	-k INT, --kmersize=INT
						kmer size (default: 31)
	-a INT, --abundance=INT
						minimum abundance (default: 1)
	-t INT, --threads=INT
    					number of threads, in [1..16]  (default: 1)
	-g STRING, --genome-type=STRING
						genome type: linear|circular (default: linear)
	-b, --build-only    build the de bruijn graph and then exit

# Removing non-maximal omnitigs

The outputted unitigs and Y-to-V contigs are maximal, but omnitigs are not necessarily. To remove non-maximal omnitigs, run for example
	
	bin/maximality chr10.k55.circular.omnitigs