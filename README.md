# complete-contigs
An assembler of non-switching and complete contigs

Usage: 
	complete-contigs OPTIONS

Brief example:
	complete-contigs -i <input file without .fa extension>  -k 31 -a 2 [-t 4] [--fastq]

Options:
	-h, --help          show this help message and exit
	--version           show program's version number and exit
	-i STRING, --input=STRING
						input file, without extension (default assumed extension .fa)
	-q, --fastq			add this flag if the input is a fastq file (then both the reads and their 
						reverse complements are added to the graph)
	-k INT, --kmersize=INT
						kmer size (default: 31)
	-a INT, --abundance=INT
						minimum abundance (default: 1)
	-t INT, --threads=INT
    					number of threads, in [1..16]  (default: 1)
	-g STRING, --genome-type=STRING
						genome type: linear|circular (default: linear)
	-b, --build-only    build the de bruijn graph and then exit
