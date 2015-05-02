#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>

#include "rlcsa/rlcsa.h"

using namespace CSA;
using namespace std;

typedef struct {
	string nodes;
	string sequence;
} contig_struct_t;


// Get current date/time, format is YYYY-MM-DD HH:mm:ss
inline string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X\n", &tstruct);

    return buf;
}

int get_data(const string& contigFileName, 
	uchar*& data, 
	vector<contig_struct_t>& contigs,
	uint64_t& char_count
	)
{
	try 
	{
		ifstream contigFile;
		contigFile.open(contigFileName);
	   	char_count = 0;
	   	string line;

	   	string genome;
	   	// counting total number of reads and their length
	   	while (getline(contigFile , line))
	   	{
	   		line = line.substr(2,line.length() - 2); // removing the leading "> "
	   		char_count += line.length() + 1; // +1 for the '\0' separator
	   		contig_struct_t contig_struct;
	   		contig_struct.nodes = line;
	   		getline(contigFile , contig_struct.sequence); // this is the actual string
	   		contigs.push_back(contig_struct);
	   	}

	   	contigFile.close();

	   	cout << "Input file " << contigFileName << " contains " << contigs.size() << " contigs." << endl;
	   	cout << "The temporary data array will have size " << char_count/(double)1000000 << "MB." << endl;

	   	uint64_t i = 0;
		data = new uchar[char_count];
		for (auto contig_struct : contigs)
		{
	    	for (uint64_t j = 0; j < contig_struct.nodes.length(); j++)
	    	{
	    		data[i] = (uchar)contig_struct.nodes[j];
	    		i++;
	    		
	    	}
	    	data[i] = '\0';
	    	i++;		
		}
	   	cout << "Created the data array" << endl;

	   	return EXIT_SUCCESS;
	} catch (exception& error)
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

}


void print_maximal_contigs(const string& contigFileName, 
	vector<contig_struct_t>& contigs,
	RLCSA& rlcsa)
{
	unordered_set<string> already_printed_contigs;

	try 
	{
		ofstream contigFile;
		contigFile.open(contigFileName + ".maximal");

		// looping over the contigs and
		// checking if they appear twice in the rlcsa or not
		for (auto contig_struct : contigs)
		{
			bool left_maximal = true;
			bool right_maximal = true;
			pair_type result_left = rlcsa.count("," + contig_struct.nodes);
			if (length(result_left) > 0)
			{
				left_maximal = false;
			}
			pair_type result_right = rlcsa.count(contig_struct.nodes + ",");
			if (length(result_right) > 0)
			{
				right_maximal = false;
			}

			if (left_maximal and right_maximal) 
			{ // maximal, but may ocurr more times
				pair_type result = rlcsa.count(contig_struct.nodes);
				if (length(result) == 1)
				{
					contigFile << "> " << contig_struct.nodes << endl;
					contigFile << contig_struct.sequence << endl;
				} 
				else if (already_printed_contigs.count(contig_struct.nodes) == 0)
				{
					already_printed_contigs.insert(contig_struct.nodes);
					contigFile << "> " << contig_struct.nodes << endl;
					contigFile << contig_struct.sequence << endl;	
				}
			}
		}
		contigFile.close();
	} 
	catch (exception& error)
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
	}
}

int main(int argc, char** argv)
{
    uint64_t char_count;
    uchar *data = NULL;
    string contigFileName;
	vector<contig_struct_t> contigs;

 	double startTime = readTimer();

	if (argc <= 1) 
	{
		cerr << "Usage: ./maximality <contig file> " << endl;
	 	return 0;
	}

	contigFileName = argv[1];

	cout << contigFileName << endl;

	if (EXIT_FAILURE == get_data(contigFileName, data, contigs, char_count))
	{
		return EXIT_FAILURE;
	}
	// Build RLCSA and report some information.
	RLCSA rlcsa(data, char_count, 32, 0, 1, true); // parameter with value '1' is number of threads; available is compiles with muti-thread support
	data = 0; // The constructor deleted the data.

	if ( !(rlcsa.isOk()) ) 
	{
		return EXIT_FAILURE;
	}
	rlcsa.printInfo();
	rlcsa.reportSize(true);
	rlcsa.writeTo(contigFileName);
	
	cout << "Time for loading the data and building the RLCSA: " << readTimer() - startTime << "sec" << endl;
 	
 	// we sample internal nodes and check their abundances
 	startTime = readTimer();

 	print_maximal_contigs(contigFileName, contigs, rlcsa);

 	cout << "Time for printing maximal contigs: " << readTimer() - startTime << "sec" << endl;
 	cout << "Time: " << currentDateTime();

	return EXIT_SUCCESS;
}