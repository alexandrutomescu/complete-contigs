
#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <omp.h>
#include <lemon/list_graph.h>
#include <lemon/static_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/maps.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <assert.h>
#include <time.h>
// #include <getopt.h>

#include "OptionParser.h"

#define GET_GRAPH(tid) (tid == 1) ? g1 : (tid == 2) ? g2 : (tid == 3) ? g3 : (tid == 4) ? g4 : (tid == 5) ? g5 : (tid == 6) ? g6 : (tid == 7) ? g7 : (tid == 8) ? g8 : (tid == 9) ? g9 : (tid == 10) ? g10 : (tid == 11) ? g11 : (tid == 12) ? g12 : (tid == 13) ? g13 : (tid == 14) ? g14 : (tid == 15) ? g15 : g16
//#define GET_NODE(tid,node) (tid == 1) ? nm1[node] : (tid == 2) ? nm2[node] : (tid == 3) ? nm3[node] : (tid == 4) ? nm4[node] : (tid == 5) ? nm5[node] : (tid == 6) ? nm6[node] : (tid == 7) ? nm7[node] : (tid == 8) ? nm8[node] : (tid == 9) ? nm9[node] : (tid == 10) ? nm10[node] : (tid == 11) ? nm11[node] : (tid == 12) ? nm12[node] : (tid == 13) ? nm13[node] : (tid == 14) ? nm14[node] : (tid == 15) ? nm15[node] : nm16[node]
#define GET_NODE_MAP(tid) (tid == 1) ? nm1 : (tid == 2) ? nm2 : (tid == 3) ? nm3 : (tid == 4) ? nm4 : (tid == 5) ? nm5 : (tid == 6) ? nm6 : (tid == 7) ? nm7 : (tid == 8) ? nm8 : (tid == 9) ? nm9 : (tid == 10) ? nm10 : (tid == 11) ? nm11 : (tid == 12) ? nm12 : (tid == 13) ? nm13 : (tid == 14) ? nm14 : (tid == 15) ? nm15 : nm16

extern int N_THREADS;

using namespace std;
using namespace lemon;

struct contig {
	std::list<int> nodes;
	std::string str;
};

struct contig_v {
	std::vector<int> nodes;
	std::string str;
};

typedef pair<int,int> pair_of_ints;

size_t hash_pair( const pair_of_ints& p );
typedef unordered_set<pair_of_ints,function<decltype(hash_pair)>> set_of_pairs;

void make_upper_case(string& s);

inline string reverse_complement_char(const char c)
{
	if (c == 'A') return "T";
	if (c == 'T') return "A";
	if (c == 'C') return "G";
	if (c == 'G') return "C";
	return "N";
}

inline string reverse_complement(const string& s)
{
	string ret = "";
	for (uint64_t i = 1; i <= s.length(); i++)
	{
		ret += reverse_complement_char(s[s.length() - i]);
	}
	return ret;
}

// Check if a file is readable
bool is_readable( const std::string & file );

size_t count_unary_nodes(StaticDigraph& graph);

size_t count_unary_arcs(StaticDigraph& graph);

size_t count_unary_arcs_ld(ListDigraph& graph);

inline bool is_unary_node(ListDigraph& graph, ListDigraph::Node& node)
{
	if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) == 1))
		return true;
	else
		return false;
}

inline string plus_strings(const string& a, const string& b, size_t kmersize)
{
	if (a == "") return b;
	if (b == "") return a;
	string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
	return ret;
}

void populate_with_strings(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const StaticDigraph::NodeMap<size_t>& length,
	vector<contig>& collection);

void populate_with_strings_from_node_labels(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection);

void populate_with_strings_list_digraph(const string& sequence, 
	size_t kmersize, 
	const ListDigraph& graph,
	const ListDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection);

void print_collection(vector<contig>& collection, string inputFileName, string suffix);

void compute_statistics(vector<contig>& collection, size_t seqLength);

void save_safe_pairs(set_of_pairs& safe_pairs, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type);

void try_loading_omnitigs(vector<contig>& omnitigs,
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type);

void print_graph_in_dot(StaticDigraph& graph, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type);

int load_data(string& sequence, 
	size_t& seqLength, 
	const string& inputFileName,
	const size_t kmersize,
	const bool circular_genome,
	const bool do_not_contract_arcs,
	StaticDigraph& graph,
	StaticDigraph::NodeMap<size_t>& length,
	StaticDigraph::NodeMap<size_t>& seqStart,
	StaticDigraph::NodeMap<string>& nodeLabel,
	set_of_pairs& safe_pairs
);

int load_data_from_reads(string& sequence, 
	size_t& seqLength, 
	const string& inputFileName,
	const size_t kmersize,
	const bool circular_genome,
	const bool do_not_contract_arcs,
	const int abundance,
	StaticDigraph& graph,
	StaticDigraph::NodeMap<size_t>& length,
	StaticDigraph::NodeMap<size_t>& seqStart,
	set_of_pairs& safe_pairs
);

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

#endif // UTILS_H_INCLUDED
