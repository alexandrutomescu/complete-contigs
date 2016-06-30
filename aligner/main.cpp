// g++ -Wall -O3 -std=c++0x -DMASSIVE_DATA_RLCSA -o aligner main.cpp rlcsa/rlcsa.a

#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <sstream>
#include <algorithm>

#include "rlcsa/rlcsa.h"
#include "OptionParser.h"

size_t GENOME_LENGTH;

using namespace CSA;
using namespace std;

struct SNP_info 
{
	string contig_name;
	int64_t pos_on_contig;
	int64_t max_snps_on_contig;
};


uint64_t get_genome_length(string referenceFileName)
{
	ifstream genomeFile;
	genomeFile.open(referenceFileName);

	string genome;
	string line;
   	getline(genomeFile , line); // this is the comment line
   	while (getline(genomeFile , line))
   	{
   		genome += line;
   	}

   	genomeFile.close();

    return genome.length();
}

void make_upper_case(string& s)
{
	for (size_t i = 0; i < s.length(); i++)
	{
		if (s[i] == 'a') s[i] = 'A';
		else if (s[i] == 'c') s[i] = 'C';
		else if (s[i] == 'g') s[i] = 'G';
		else if (s[i] == 't') s[i] = 'T';
		else if (s[i] == 'n') s[i] = 'N';
		else if (s[i] == 'A') s[i] = 'A';
		else if (s[i] == 'C') s[i] = 'C';
		else if (s[i] == 'G') s[i] = 'G';
		else if (s[i] == 'T') s[i] = 'T';
		else s[i] = 'N';
	}
}

string reverse_complement_char(const char c)
{
	if (c == 'A') return "T";
	if (c == 'T') return "A";
	if (c == 'C') return "G";
	if (c == 'G') return "C";
	if (c == 'N') return "N";

	return "N";
	// if (c.compare("A") == 0) return "T";
	// if (c.compare("T") == 0) return "A";
	// if (c.compare("C") == 0) return "G";
	// if (c.compare("G") == 0) return "C";
	// if (c.compare("N") == 0) return "N";
}

string reverse_complement(const string& s)
{
	string ret = "";
	for (uint i = 1; i <= s.length(); i++)
	{
		ret += reverse_complement_char(s[s.length() - i]);
	}
	return ret;
}

int get_reads(const string readFileName, 
	vector<string>& reads,
	vector<string>& read_labels
	)
{
	try 
	{
		ifstream readFile;
		readFile.open(readFileName);
	   	string line;

	   	while (getline(readFile , line)) // this is the comment line
	   	{
	   		read_labels.push_back(line);
	    	getline(readFile , line); // the actual read
	    	reads.push_back(line);
	   	}
	   	readFile.close();
	} catch (exception& error) 
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

	cout << "Input file contains " << reads.size() << " reads." << endl;

   	return EXIT_SUCCESS;
}

int get_data_for_rlcsa(const string& genomeFileName, 
	bool circular_genome,
	uchar*& data, 
	uint64_t& char_count
	)
{
	ifstream readFile;
	readFile.open(genomeFileName);
   	char_count = 0;
   	string line;
   	vector<string> reads;

   	string genome;
   	// counting total number of reads and their length
   	getline(readFile , line); // this is the comment line
   	while (getline(readFile , line))
   	{
   		genome += line;
   	}
    make_upper_case(genome);

    GENOME_LENGTH = genome.length();

    if (circular_genome)
    {
    	genome = genome + genome;
    }

    reads.push_back(genome);
    char_count = char_count + (genome.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'

    // /// REVERSE COMPLEMENT
    // string genome_RC = reverse_complement(genome);
    // reads.push_back(genome_RC);
    // char_count = char_count + (genome_RC.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
    // /// REVERSE COMPLEMENT

  	reads.push_back("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
	char_count = char_count + (60 + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
   	readFile.close();

   	cout << "Input file " << genomeFileName << " contains " << reads.size() << " reads." << endl;
   	cout << "The temporary data array will have size " << char_count/(double)1000000000 << "GB." << endl;

   	uint64_t i = 0;
	data = new uchar[char_count];
	for (auto read : reads)
	{
    	for (uint64_t j = 0; j < read.length(); j++)
    	{
    		data[i] = (uchar)read[j];
    		i++;
    		
    	}
    	data[i] = '\0';
    	i++;		
	}
   	cout << "Created the data array" << endl;

   	return EXIT_SUCCESS;
}

void load_snp_info(string dbsnpFileName, 
	string chr, 
	vector<size_t>& prefixSum_snp,
	vector<size_t>& next_snp,
	unordered_map<uint64_t,SNP_info>& SNP_map,
	string referenceFileName)
{
	GENOME_LENGTH = get_genome_length(referenceFileName) + 1000000;

	vector<size_t> indicator_snp(GENOME_LENGTH,0);

	cout << "Got here" << endl;

	ifstream dbsnpFile;
	dbsnpFile.open(dbsnpFileName);
	string line;
	size_t n_SNPs = 0;


	while (getline(dbsnpFile, line))
	{
		stringstream linestream(line);
		string CHR, REF, ALT, ID;
		size_t POS;

		getline(linestream, CHR, '\t'); // CHR field
		linestream >> POS; // POS field
		getline(linestream, ID, '\t'); // '\t' after POS
		getline(linestream, ID, '\t'); // ID field
		getline(linestream, REF, '\t');
		getline(linestream, ALT, '\t');


		if ((REF.length() == 1) and (ALT.length() == 1))
		{
			// cout << CHR << ":" << POS << ":" << ID << ":" << REF << ":" << ALT << endl;
			assert(POS < GENOME_LENGTH);
			indicator_snp[POS] = 1;
			SNP_info new_snp_info;
			new_snp_info.contig_name = "";
			new_snp_info.pos_on_contig = -1;
			new_snp_info.max_snps_on_contig = -1;

			SNP_map[POS] = new_snp_info;
			n_SNPs++;
			// cout << "Setting a SNP on position " << POS << " as " << indicator_snp[POS] << endl;
		}
	}

	cout << "Read all SNPs " << endl;

	prefixSum_snp.resize(GENOME_LENGTH,0);
	prefixSum_snp[0] = 0;

	for (size_t i = 1; i < GENOME_LENGTH; i++)
	{
		prefixSum_snp[i] = indicator_snp[i - 1] + prefixSum_snp[i - 1];
		if (prefixSum_snp[i] > prefixSum_snp[i - 1])
		{
			// cout << "prefixSum_snp[" << i << "]=" << prefixSum_snp[i] << endl;
		}
	}

	cout << "Computed prefixSum_snp " << endl;

	next_snp.resize(GENOME_LENGTH,0);
	next_snp[GENOME_LENGTH - 1] = (indicator_snp[GENOME_LENGTH - 1]) ? GENOME_LENGTH - 1 : GENOME_LENGTH;
	cout << "next_snp[GENOME_LENGTH - 1] = " << next_snp[GENOME_LENGTH - 1] << endl;

	cout << "Starting the for loop for next_snp " << endl;
	for (int64_t i = GENOME_LENGTH - 2; i >= 0; i--)
	{
		if (indicator_snp[i] == 1) 
		{
			next_snp[i] = i;
		}
		else
		{
			next_snp[i] = next_snp[i + 1];
		}
		// cout << "next_snp[" << i << "] = " << next_snp[i] << endl;
	}

	cout << "GENOME_LENGTH " << GENOME_LENGTH << endl;
	cout << "Found n_SNPs = " << n_SNPs << endl;

	dbsnpFile.close();
}

void load_gtf_info(string gtfFileName, 
	vector<size_t>& next_exon_start,
	vector<size_t>& next_exon_end,
	string referenceFileName)
{
	GENOME_LENGTH = get_genome_length(referenceFileName) + 1000000;

	next_exon_start.resize(GENOME_LENGTH,0);
	next_exon_end.resize(GENOME_LENGTH,0);
	vector<size_t> indicator_start(GENOME_LENGTH,0);
	vector<size_t> indicator_end(GENOME_LENGTH,0);

	ifstream gtfFile;
	gtfFile.open(gtfFileName);
	string line;

	while (getline(gtfFile, line))
	{
		stringstream linestream(line);
		string CHR, NAME, EXON;
		size_t START, END;

		getline(linestream, CHR, '\t'); 
		getline(linestream, NAME, '\t'); 
		getline(linestream, EXON, '\t');

		linestream >> START;
		linestream >> END;

		assert(END < GENOME_LENGTH);
		indicator_start[START] = 1;
		indicator_end[END] = 1;
	}

	cout << "Read all exons " << endl;

	next_exon_start[GENOME_LENGTH - 1] = (indicator_start[GENOME_LENGTH - 1]) ? GENOME_LENGTH - 1 : GENOME_LENGTH;
	next_exon_end[GENOME_LENGTH - 1] = (indicator_end[GENOME_LENGTH - 1]) ? GENOME_LENGTH - 1 : GENOME_LENGTH;

	for (int64_t i = GENOME_LENGTH - 2; i >= 0; i--)
	{
		if (indicator_start[i] == 1) 
		{
			next_exon_start[i] = i;
		}
		else
		{
			next_exon_start[i] = next_exon_start[i + 1];
		}
		////////////////////////////
		if (indicator_end[i] == 1) 
		{
			next_exon_end[i] = i;
		}
		else
		{
			next_exon_end[i] = next_exon_end[i + 1];
		}
	}

	gtfFile.close();
}


void align_reads(const RLCSA* rlcsa,
	vector<string>& reads,
	vector<string>& read_labels,
	ofstream& fileStats,
	ofstream& fileDumpContigStats,
	ofstream& fileDumpSNPStats,
	ofstream& fileReads,
	string referenceFileName,
	string previous_line,
	vector<size_t>& prefixSum_snp,
	vector<size_t>& next_snp,
	unordered_map<uint64_t,SNP_info> SNP_map,
	vector<size_t>& next_exon_start,
	vector<size_t>& next_exon_end)
{

	bool snp_analysis = prefixSum_snp.size();
	bool exon_analysis = next_exon_end.size() > 0 ? true : false;

	fileStats << previous_line << ",";

	uint64_t genomeLength = get_genome_length(referenceFileName);
	double totalLength = 0;
	uint64_t total_snps = 0;
	uint64_t total_exon_content = 0;
	uint64_t total_n_exons = 0;
	double maxLength = 0;
	int64_t max_snps_per_contig = 0;
	uint64_t notMapping = 0;
	uint64_t n_aligned_reads = 0;

	double n_SNP_pairs = 0;

	cout << "genomeLength: " << genomeLength << endl;
	fileStats << genomeLength << ",";

	// vector<int> readsStartingAtPosition[genomeLength];
	// vector<int> readsEndingAtPosition[genomeLength];
	vector< vector<int> > readsStartingAtPosition(genomeLength, vector<int>());
	vector< vector<int> > readsEndingAtPosition(genomeLength, vector<int>());

	cout << "Starting to align the contigs" << endl;

	for (usint i = 0; i < reads.size(); i++)
	{
		fileReads << read_labels[i] << endl;
		fileReads << reads[i];
		pair_type result = rlcsa->count(reads[i]);
		if (length(result) > 0)
		{
			n_aligned_reads++;
			int64_t max_SNPs_in_alignment = 0;
			uint64_t max_exon_content = 0;
			uint64_t max_n_exons = 0;
			// updating the stats about ContigNameAndPos
			string contig_name = read_labels[i];
			replace( contig_name.begin(), contig_name.end(), ',', '_');
			contig_name = contig_name.substr(2);

			fileReads << " maps" << endl;
			usint* positions = rlcsa->locate(result);

			// check that it is not all made up of 'N's
			bool only_Ns = true;
			for (uint64_t j = 0; j < reads[i].length(); j++)
			{
				if (reads[i][j] != 'N')
				{
					only_Ns = false;
					break;
				}
			}

			for (usint j = 0; (not only_Ns) and (j < length(result)); j++)
			{
				if (positions[j] < genomeLength) // for circular genomes
				{
					readsStartingAtPosition[positions[j]].push_back(i);	

					uint64_t x = positions[j];
					uint64_t y = positions[j] + reads[i].length();

					if (snp_analysis)
					{
						int64_t SNPs_in_alignment = prefixSum_snp[y] - prefixSum_snp[x];
						if (max_SNPs_in_alignment < SNPs_in_alignment)
						{
							max_SNPs_in_alignment = SNPs_in_alignment;
						}

						uint64_t z = next_snp[x];
						// z scans the positions of the alignment where there are SNPs
						while (z < y)
						{	
							if (SNP_map[z].max_snps_on_contig < SNPs_in_alignment)
							{
								SNP_map[z].max_snps_on_contig = SNPs_in_alignment;
								SNP_map[z].contig_name = contig_name;
								SNP_map[z].pos_on_contig = z - x;
							}
							z = next_snp[z + 1];
						}
					}

					if (exon_analysis)
					{
						uint64_t exon_content = 0;
						uint64_t n_exons = 0;
						uint64_t z_start = next_exon_start[x];
						uint64_t z_end = next_exon_end[x];

						// if exon has started already
						if (z_end < z_start) 
						{
							exon_content += z_end - x + 1;
							n_exons++;
							z_end = next_exon_end[z_start];
						}

						while ((z_start < y) and (z_end < y))
						{
							exon_content += z_end - z_start + 1;
							n_exons++;

							z_start = next_exon_start[z_end + 1];
							z_end = next_exon_end[z_start];
						}

						// if exon ends after contig end
						if (z_start < y)
						{
							exon_content += y - z_start + 1;
							n_exons++;
						}

						if (exon_content > max_exon_content)
						{
							max_exon_content = exon_content;
							max_n_exons = n_exons;
						}
					}
				}
				if (positions[j] + reads[i].length() < genomeLength)
				{
					readsEndingAtPosition[positions[j] + reads[i].length()].push_back(i);
				}
			}


			if (snp_analysis)
			{
				// one line per contig: contig length, max #SNPs
				total_snps += max_SNPs_in_alignment;
				n_SNP_pairs += max_SNPs_in_alignment * (max_SNPs_in_alignment - 1) / 2;
				if (max_SNPs_in_alignment > max_snps_per_contig)
				{
					max_snps_per_contig = max_SNPs_in_alignment;
				}
			}
			if (exon_analysis)
			{
				total_exon_content += max_exon_content;
				total_n_exons += max_n_exons;
			}

			fileDumpContigStats << contig_name << "," << reads[i].length() << "," << max_SNPs_in_alignment << "," << max_n_exons << "," << max_exon_content << "," << length(result) << endl;

			// totalLength and maxLength
			totalLength += reads[i].length();
			if (reads[i].length() > maxLength)
			{
				maxLength = reads[i].length();
			}
		}
		else
		{
			fileReads << " -1" << endl;
			notMapping++;
		}
	}

	fileReads.close();

	multiset<int> active_reads;
	double sum = 0;
	for (usint i = 0; i < genomeLength; i++)
	{
		for (auto read_idx : readsStartingAtPosition[i])
		{
			active_reads.insert(read_idx);
		}
		for (auto read_idx : readsEndingAtPosition[i])
		{
			active_reads.erase(read_idx);
		}
		double curr_sum = 0;

		if (active_reads.size() > 0)
		{
			for (auto read_idx : active_reads)
			{
				curr_sum = curr_sum + reads[read_idx].length();
			}

			double curr_avg = curr_sum / active_reads.size();
			sum = sum + curr_avg;
		}

	}

	double esize = sum / genomeLength;


	if (snp_analysis)
	{
		// one line per SNP X: pos, contig name covering X and covering most SNPs, pos of X on contig, #SNPs in the same contig
		uint64_t z = next_snp[0];
		while (z < GENOME_LENGTH)
		{
			fileDumpSNPStats << z << "," << SNP_map[z].contig_name << "," << SNP_map[z].pos_on_contig << "," << SNP_map[z].max_snps_on_contig << endl;
			z = next_snp[z + 1];
		}
		
	}


	// outputStatsFile << "#Contigs,AVG length,MAX length,E-size,Unmapped" << endl;
	
	fileStats << reads.size() << "," << totalLength / n_aligned_reads << "," << maxLength << "," << esize << "," << (float)total_snps / n_aligned_reads << "," << max_snps_per_contig << "," << n_SNP_pairs << "," << total_exon_content << "," << (float)total_exon_content / n_aligned_reads << "," << total_n_exons << "," << (float)total_n_exons / n_aligned_reads << "," << notMapping << endl;

	cout << "********* STATS *********** " << endl;
	cout << "#Contigs   : " << reads.size() << endl;
	cout << "Unmapped   : " << notMapping << endl;
	cout << "AVG length : " << totalLength / n_aligned_reads << endl;
	cout << "MAX length : " << maxLength << endl;
	cout << "E-size     : " << esize << endl;
	cout << "AVG #SNPs/c: " << total_snps / n_aligned_reads << endl;
	cout << "MAX #SNPs/c: " << max_snps_per_contig << endl;
	cout << "total exon content  : " << total_exon_content << endl;
	cout << "AVG exon content    : " << (float)total_exon_content / n_aligned_reads << endl;
	cout << "total n_exons  	 : " << total_n_exons << endl;
	cout << "AVG n_exons         : " << (float)total_n_exons / n_aligned_reads << endl;
	cout << "*************************** " << endl;

}

int main(int argc, char** argv)
{
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName;
    string referenceFileName;
    string outputStatsFileName;
    string dbsnpFileName = "";
    string gtfFileName = "";
    string chr;
    bool buildindex = false;
    bool loadindex = false;
	vector<string> reads;
	vector<string> read_labels;
	bool circular_genome = false;

	// command line argument parser
	string usage = "\n  %prog OPTIONS"
		"\n\nBrief example:"
		"\n  Usage: ./aligner [--loadindex|--buildindex] [--circular] --reference <.fa file> --readfile <.fastq file> --outputfile <.csv file> [--dbsnp dbsnp.vcf] [--chr 14]";
	const string version = "%prog 0.2\nCopyright (C) 2014-2015 Alexandru Tomescu\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>.\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.";
	const string desc = "";
	const string epilog = "";
	
	optparse::OptionParser parser = optparse::OptionParser()
    	.usage(usage)
    	.version(version)
    	.description(desc)
    	.epilog(epilog);

	parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("contig file (.fasta/.fa format)");
	parser.add_option("-g", "--reference") .type("string") .dest("g") .set_default("") .help("reference genome (default: %default)");
	parser.add_option("-b", "--buildindex") .action("store_true") .set_default(false) .dest("b") .help("build the index and exit");
	parser.add_option("-l", "--loadindex") .action("store_true") .set_default(false) .dest("l") .help("load the index from disk");
	parser.add_option("-c", "--circular") .action("store_true") .set_default(false) .dest("c") .help("flag for circular genome");
	parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output csv file [the one generated by the assembler]");
	parser.add_option("-d", "--dbsnp") .type("string") .dest("d") .set_default("") .help("dbsnp VCF file");
	parser.add_option("-e", "--gtf") .type("string") .dest("e") .set_default("") .help("gene annotation GTF file");
	// parser.add_option("-x", "--chr") .type("string") .dest("x") .set_default("") .help("chromosome number");

	optparse::Values& options = parser.parse_args(argc, argv);

	readFileName 		= (string) options.get("r");
	referenceFileName 	= (string) options.get("g");
	buildindex 			= (options.get("b") ? true : false);
	loadindex 			= (options.get("l") ? true : false);
	circular_genome		= (options.get("c") ? true : false);
	outputStatsFileName = (string) options.get("o");
	dbsnpFileName 		= (string) options.get("d");
	gtfFileName 		= (string) options.get("e");
	// chr 				= (string) options.get("x");

	if (((buildindex == false) and (loadindex == false)) or ((buildindex == true) and (loadindex == true)))
	{
		cerr << "You must specify either --buildindex or --loadindex, and not both" << endl;
		return EXIT_FAILURE;
	}

	cout << referenceFileName << endl;
	cout << "circular_genome: " << circular_genome << endl;

	// we need to build the index and exit
	if (buildindex) 
	{
		if (EXIT_FAILURE == get_data_for_rlcsa(referenceFileName, circular_genome, data, char_count))
		{
			return EXIT_FAILURE;
		}
		// Build RLCSA and report some information.
		RLCSA rlcsa(data, char_count, 32, 32, 1, true); // parameter with value '1' is number of threads; available if compiled with muti-thread support
 		data = 0; // The constructor deleted the data.

		if ( !(rlcsa.isOk()) ) 
 		{
 			return EXIT_FAILURE;
 		}
 		rlcsa.printInfo();
 		rlcsa.reportSize(true);
 		rlcsa.writeTo(referenceFileName);
 		return EXIT_SUCCESS;
	}

	// we need to load the index
 	const RLCSA* rlcsa = new RLCSA(referenceFileName, false);
 	if (!(rlcsa->isOk())) 
 	{
 		return EXIT_FAILURE;
 	}
 	rlcsa->printInfo();
 	rlcsa->reportSize(true);

 	cout << "Start..." << endl;

 	vector<size_t> prefixSum_snp; // prefixSum_snp[i] := number of SNPs from 0..i-1
 	vector<size_t> next_snp; // prefixSum_snp[i] := number of SNPs from 0..i-1
 	unordered_map<uint64_t,SNP_info> SNP_map;
	vector<size_t> next_exon_start;
	vector<size_t> next_exon_end;

 	if (dbsnpFileName != "")
 	{
 		cout << "Loading prefix sums" << endl;
 		load_snp_info(dbsnpFileName, chr, prefixSum_snp, next_snp, SNP_map, referenceFileName);
 	}

 	if (gtfFileName != "")
 	{
 		cout << "Loading GTF annotation" << endl;
 		load_gtf_info(gtfFileName, next_exon_start, next_exon_end, referenceFileName);
 	}

	ifstream fileStats_read;
	string header_line, unitigs_line, ns_contigs_line, YtoV_contigs_line, omnitigs_line, previous_line;

	fileStats_read.open(outputStatsFileName);
	getline(fileStats_read,header_line);
	getline(fileStats_read,unitigs_line);
	//getline(fileStats_read,ns_contigs_line);
	getline(fileStats_read,YtoV_contigs_line);
	getline(fileStats_read,omnitigs_line);
	fileStats_read.close();
	
	ofstream fileStats;
	fileStats.open(outputStatsFileName);
	fileStats << header_line << ",genome length,#strings,AVG_length,MAX_length,E-size,AVG #SNPs/c,MAX #SNPs/c,#SNP pairs in same contig,total exon content,AVG exon content,total n_exons,AVG n_exons,unmapped" << endl;

	for (string algorithm : { "unitigs",  "YtoV-contigs", "omnitigs.maximal"}) // "non-switching-contigs",
	{

		ofstream fileDumpContigStats;
		ofstream fileDumpSNPStats;
		fileDumpContigStats.open(readFileName + "." + algorithm + ".dumpContigStats.csv");
		fileDumpContigStats << "contig_name,contig_length,#SNPs,#exons touched,exon content,number of alignments (max over all its alignments)" << endl;

		fileDumpSNPStats.open(readFileName + "." + algorithm + ".dumpSNPStats.csv");
		fileDumpSNPStats << "pos,contig_name,pos_on_contig,#SNPs in the same contig (max over all contigs aligning on that SNP)" << endl;

		cout << "*********** " << algorithm << endl;

		if (algorithm == "unitigs")
		{
			previous_line = unitigs_line;
		}
		if (algorithm == "non-switching-contigs")
		{
			previous_line = ns_contigs_line;
		}
		if (algorithm == "YtoV-contigs")
		{
			previous_line = YtoV_contigs_line;
		}		
		if (algorithm == "omnitigs.maximal")
		{
			previous_line = omnitigs_line;
		}

	 	// we load the reads
 		get_reads(readFileName + "." + algorithm, reads, read_labels);
 		
 		ofstream fileReads;	
		fileReads.open(readFileName + "." + algorithm + ".aligned");
 		align_reads(rlcsa,
			reads,
			read_labels,
			fileStats,
			fileDumpContigStats,
			fileDumpSNPStats,
			fileReads,
			referenceFileName,
			previous_line,
			prefixSum_snp,
			next_snp,
			SNP_map,
			next_exon_start,
			next_exon_end);

 		reads.clear();
 		read_labels.clear();

 		fileDumpContigStats.close();
 		fileDumpSNPStats.close();
 	}

 	fileStats.close();


	return EXIT_SUCCESS;
}