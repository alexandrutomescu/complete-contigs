// compile command:
// on my mac: g++ main.cpp -o ec-contigs -std=c++0x -I ./lemon-1.3.1/include -O3
// on linux: g++ main.cpp -o ec-contigs -std=c++0x -I ~/local/include/ -O3 -fopenmp

#include "utils.h"
//#include "maximality.h"

using namespace std;
using namespace lemon;

int N_THREADS;


vector<contig> compute_unitigs(const StaticDigraph& graph, 
	const StaticDigraph::NodeMap<size_t>& length, 
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;
	unordered_set<int> processedArcs;

	// first, merging the arcs incident to unary nodes
	for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) == 1))
		{
			StaticDigraph::InArcIt in_arc(graph,node);
			StaticDigraph::OutArcIt out_arc(graph,node);

			contig entry;
			entry.nodes.push_back(graph.id(graph.source(in_arc)));
			entry.nodes.push_back(graph.id(node));
			entry.nodes.push_back(graph.id(graph.target(out_arc)));

			ret.push_back(entry);

			processedArcs.insert(graph.id(in_arc));
			processedArcs.insert(graph.id(out_arc));
		}
	}

	for (StaticDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (processedArcs.count(graph.id(a)) == 0)	
		{
			
			contig entry;
			entry.nodes.push_back(graph.id(graph.source(a)));
			entry.nodes.push_back(graph.id(graph.target(a)));
			ret.push_back(entry);	

		}
	}

	return ret;
}

vector<contig> compute_nice_contigs(const StaticDigraph& graph, 
	const StaticDigraph::NodeMap<size_t>& length, 
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;

	cout << "Entered nice_contigs" << endl;

	StaticDigraph::Node current_node;
	string s;

	// iterating over all arcs
	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)	
	{
		if (((countOutArcs(graph,graph.source(arc)) > 1) and (countInArcs(graph,graph.target(arc)) > 1)) or 
		 	((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.source(arc)) == 1)))
		{
			// traversing the outtig (forward)
			contig entry;

			//cout << "creating a nice contig form arc " << graph.id(graph.source(arc)) << "," << graph.id(graph.target(arc)) << endl;
			current_node = graph.target(arc);
			while (countOutArcs(graph, current_node) == 1)
			{
				entry.nodes.push_back(graph.id(current_node));
				StaticDigraph::OutArcIt out_arc(graph,current_node);
				if (current_node == graph.target(out_arc))
				{
					entry.nodes.push_back(graph.id(current_node));
					break;
				}
				current_node = graph.target(out_arc);
				//cout << "seen " << graph.id(current_node) << " forward" << endl;
			}
			entry.nodes.push_back(graph.id(current_node));

			// traversing the intig (backward)
			current_node = graph.source(arc);
			while (countInArcs(graph, current_node) == 1)
			{
				entry.nodes.push_front(graph.id(current_node));
				StaticDigraph::InArcIt in_arc(graph,current_node);
				if (current_node == graph.source(in_arc))
				{
					entry.nodes.push_front(graph.id(current_node));
					break;
				}
				current_node = graph.source(in_arc);
				//cout << "seen " << graph.id(current_node) << " backward" << endl;
			}
			entry.nodes.push_front(graph.id(current_node));

			ret.push_back(entry);
		}		
	}

	return ret;
}

vector<contig> compute_ec_contigs_2(StaticDigraph& graph, 
	StaticDigraph::NodeMap<size_t>& length, 
	StaticDigraph::NodeMap<size_t>& seqStart,
	set_of_pairs& safe_pairs,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;
	size_t n_nodes = countNodes(graph);
	ret.reserve(n_nodes * 100);

	cout << "Entered compute_ec_contigs_2" << endl;


	StaticDigraph g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16;
	StaticDigraph::NodeMap<StaticDigraph::Node> nm1(graph),nm2(graph),nm3(graph),nm4(graph),nm5(graph),nm6(graph),nm7(graph),nm8(graph),nm9(graph),nm10(graph),nm11(graph),nm12(graph),nm13(graph),nm14(graph),nm15(graph),nm16(graph);


	switch (N_THREADS) {
		case 16: digraphCopy(graph,g16).nodeRef(nm16).run();
		case 15: digraphCopy(graph,g15).nodeRef(nm15).run();
		case 14: digraphCopy(graph,g14).nodeRef(nm14).run();
		case 13: digraphCopy(graph,g13).nodeRef(nm13).run();
		case 12: digraphCopy(graph,g12).nodeRef(nm12).run();
		case 11: digraphCopy(graph,g11).nodeRef(nm11).run();
		case 10: digraphCopy(graph,g10).nodeRef(nm10).run();
		case 9: digraphCopy(graph,g9).nodeRef(nm9).run();
		case 8: digraphCopy(graph,g8).nodeRef(nm8).run();
		case 7: digraphCopy(graph,g7).nodeRef(nm7).run();
		case 6: digraphCopy(graph,g6).nodeRef(nm6).run();
		case 5: digraphCopy(graph,g5).nodeRef(nm5).run();
		case 4: digraphCopy(graph,g4).nodeRef(nm4).run();
		case 3: digraphCopy(graph,g3).nodeRef(nm3).run();
		case 2: digraphCopy(graph,g2).nodeRef(nm2).run();
		case 1: digraphCopy(graph,g1).nodeRef(nm1).run();
	}

	cout << "Finished making copies" << endl;

	size_t progress = 0;

	// for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	// #pragma omp parallel for
	omp_set_dynamic(0);
	#pragma omp parallel for num_threads(N_THREADS) shared(ret,safe_pairs)
	for (uint64_t i = 0; i < n_nodes; i++)
	{
		StaticDigraph::Node node = graph.node(i);

		//************** progress ************** //
		#pragma omp critical
		{
			progress = progress + 1;
		}

		if (progress % 1000 == 0)
		{
			cout << "Node #" << progress << "/" << n_nodes << " ";
		  	cout << "Time: " << currentDateTime();
		  	cout.flush();
		}
		//************** progress ************** //

		// if one-in, then all right extensions are safe
		if (countInArcs(graph,node) == 1)
		{
			StaticDigraph::InArcIt in_arc(graph,node);
			for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
			{
				if (out_arc != in_arc)
				{
					contig entry;
					entry.nodes.push_back(graph.id(graph.source(in_arc)));
					entry.nodes.push_back(graph.id(node));
					entry.nodes.push_back(graph.id(graph.target(out_arc)));

					#pragma omp critical
					{
						ret.push_back(entry);
						safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));	
					}	
				}
			}	
		} else 
		// if one-out, then all left extensions are safe
		if (countOutArcs(graph,node) == 1)
		{
			StaticDigraph::OutArcIt out_arc(graph,node);
			for (StaticDigraph::InArcIt in_arc(graph,node); in_arc != INVALID; ++in_arc)
			{
				if (in_arc != out_arc)
				{
					contig entry;
					entry.nodes.push_back(graph.id(graph.source(in_arc)));
					entry.nodes.push_back(graph.id(node));
					entry.nodes.push_back(graph.id(graph.target(out_arc)));

					#pragma omp critical
					{
						ret.push_back(entry);
						safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));	
					}
				}				
			}				
		} else
		{
			for (StaticDigraph::InArcIt in_arc(graph,node); in_arc != INVALID; ++in_arc)
			{
				for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
				{
					if (out_arc == in_arc)
					{
						continue;
					}
					// checking whether we the walk in_arc,out_arc forms an edge centric contig
					int tid = omp_get_thread_num();
					tid++;
					//cout << tid << endl;
					StaticDigraph::NodeMap<bool> filterNodesMap(GET_GRAPH(tid), true);
					FilterNodes<StaticDigraph> subgraph(GET_GRAPH(tid), filterNodesMap);

					// StaticDigraph::NodeMap<bool> filterNodesMap(graph, true);
					// FilterNodes<StaticDigraph> subgraph(graph, filterNodesMap);
					
					// check whether some in_nbr is reachable in G - node from some out-neighbor of node different than target(out_arc)
					subgraph.disable((GET_NODE_MAP(tid))[node]);
					Bfs<FilterNodes<StaticDigraph>> visit(subgraph);
					visit.init();

					// adding sources of the visit
					for (StaticDigraph::OutArcIt out_arc2(graph,node); out_arc2 != INVALID; ++out_arc2)
					{
						if (graph.target(out_arc2) != graph.target(out_arc))
						{
							visit.addSource((GET_NODE_MAP(tid))[graph.target(out_arc2)]);
						}
					}
					assert(visit.queueSize() > 0);

					vector<StaticDigraph::Node> bad_guys;
					// adding the other in_nbrs as the other targets of visit
					FilterNodes<StaticDigraph>::NodeMap<bool> visitTargets(subgraph,false);
					// size_t n_bad = 0;
					for (StaticDigraph::InArcIt in_arc2(graph,node); in_arc2 != INVALID; ++in_arc2)
					{
						if (in_arc2 != in_arc)
						{
							visitTargets[(GET_NODE_MAP(tid))[graph.source(in_arc2)]] = true;
							bad_guys.push_back((GET_NODE_MAP(tid))[graph.source(in_arc2)]);
						}
					}

					// if no target vertex is reached, add these arcs as ec_contig_2 and safe pair
					if (visit.start(visitTargets) == INVALID)
					{
						// check that indeed all bad guys are not reached
						// LEMON seems to have a bug and the above check is not sufficient
						bool no_bad_inbr_reached = true;
						for (size_t i = 0; i < bad_guys.size(); i++)
						{
							if (visit.reached(bad_guys[i]))
							{
								no_bad_inbr_reached = false;
								break;
							}
						}
						if (no_bad_inbr_reached)
						{
							contig entry;
							entry.nodes.push_back(graph.id(graph.source(in_arc)));
							entry.nodes.push_back(graph.id(node));
							entry.nodes.push_back(graph.id(graph.target(out_arc)));

							#pragma omp critical
							{
								ret.push_back(entry);
								safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));
							}
							
						}
					}
				}
			}
			
		}
	}

	return ret;
}


void extend_pair(StaticDigraph& graph,
	StaticDigraph& thread_graph,
	StaticDigraph::NodeMap<StaticDigraph::Node>& originalNode2ThreadNode,
	const set_of_pairs& safe_pairs,
	StaticDigraph::Arc in_arc,
	list<int> current_walk,
	vector<int> bad_in_nbrs,
	vector<contig>& collection_of_contigs,
	unordered_set<int>& reported_arcs
	)
{
	//cout << "Entered extend_pair " << endl;

	StaticDigraph::Node node = graph.target(in_arc);
	// iterate over all out-neighbors
	bool is_maximal = true;
	for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
	{
		// if the pair (in_arc,out_arc) is safe
		if (safe_pairs.count(pair_of_ints(graph.id(in_arc),graph.id(out_arc))) > 0)
		{
			// check whether from 'node', which becomes an internal vertex, there is no path back to 
			// the in-nbrs of the internal vertices of current_walk
			StaticDigraph::NodeMap<bool> filterNodesMap(thread_graph,true);
			FilterNodes<StaticDigraph> subgraph(thread_graph, filterNodesMap);
			Bfs<FilterNodes<StaticDigraph>> visit(subgraph);
			visit.init();

			// adding the sources
			for (StaticDigraph::OutArcIt out_arc2(graph,node); out_arc2 != INVALID; ++out_arc2)
			{
				if (out_arc2 != out_arc)
				{
					visit.addSource(originalNode2ThreadNode[graph.target(out_arc2)]);
				}
			}

			// adding the targets of the search as the bad_in_nbrs
			FilterNodes<StaticDigraph>::NodeMap<bool> visitTargets(subgraph,false);
			// size_t n_added = 0;
			vector<StaticDigraph::Node> bad_guys;
			for (size_t i = 0; i < bad_in_nbrs.size(); i++)
			{
				visitTargets[originalNode2ThreadNode[graph.node(bad_in_nbrs[i])]] = true;
				bad_guys.push_back(originalNode2ThreadNode[graph.node(bad_in_nbrs[i])]);
			}

			subgraph.disable(originalNode2ThreadNode[node]);
			// if no bad_in_nbr is reached
			if ((visit.queueSize() == 0) or (visit.start(visitTargets) == INVALID))
			{
				// check that indeed all bad guys are not reached
				// LEMON seems to have a bug and the above check is not sufficient
				bool no_bad_inbr_reached = true;
				for (size_t i = 0; i < bad_guys.size(); i++)
				{
					if (visit.reached(bad_guys[i]))
					{
						no_bad_inbr_reached = false;
						break;
					}
				}

				// if yes, then add the target of out_arc to current_walk
				if (no_bad_inbr_reached)
				{
					list<int> new_walk = current_walk;
					new_walk.push_back(graph.id(graph.target(out_arc)));

					// add the in-nbrs of node not on current_walk to bad_in_nbrs
					vector<int> new_bad_in_nbrs = bad_in_nbrs;
					for (StaticDigraph::InArcIt in_arc2(graph,node); in_arc2 != INVALID; ++in_arc2)
					{
						if (in_arc2 != in_arc)
						{
							new_bad_in_nbrs.push_back(graph.id(graph.source(in_arc2)));
						}
					}
					// and recurse
					is_maximal = false;
					#pragma omp critical
					{
						reported_arcs.insert(graph.id(out_arc));
					}

					// if in_arc and out_arc are arcs between the same two vertices, but in opposite directions
					// then we must stop the walk here
					// because there is no other way to stop it
					if (graph.source(in_arc) == graph.target(out_arc))
					{
						contig entry;
						entry.nodes = new_walk;
						#pragma omp critical
						{
							collection_of_contigs.push_back(entry);
						}

					} 
					else
					{
						extend_pair(graph, thread_graph, originalNode2ThreadNode, safe_pairs, out_arc, new_walk, new_bad_in_nbrs, collection_of_contigs, reported_arcs);		
					}
				}
			}
		}
	}

	// if we didn't manage to extend the contig, then it is right-maximal
	// and we add it to the collection of contigs
	// WARNING: it may not be left-maximal
	if (is_maximal)
	{
		contig entry;
		entry.nodes = current_walk;
		#pragma omp critical
		{
			collection_of_contigs.push_back(entry);
		}
	}

}


vector<contig> compute_ec_contigs(StaticDigraph& graph, 
	StaticDigraph::NodeMap<size_t>& length, 
	StaticDigraph::NodeMap<size_t>& seqStart,
	const set_of_pairs& safe_pairs
	)
{
	cout << "Entered ec_contigs" << endl;
	vector<contig> ret;
	unordered_set<int> reported_arcs;
	size_t progress = 0;

	vector<pair_of_ints> safe_pairs_v(safe_pairs.begin(),safe_pairs.end());
	int n_safe_pairs = safe_pairs_v.size();

	ret.reserve(n_safe_pairs * 100);

	StaticDigraph g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16;
	StaticDigraph::NodeMap<StaticDigraph::Node> nm1(graph),nm2(graph),nm3(graph),nm4(graph),nm5(graph),nm6(graph),nm7(graph),nm8(graph),nm9(graph),nm10(graph),nm11(graph),nm12(graph),nm13(graph),nm14(graph),nm15(graph),nm16(graph);

	switch (N_THREADS) {
		case 16: digraphCopy(graph,g16).nodeRef(nm16).run();
		case 15: digraphCopy(graph,g15).nodeRef(nm15).run();
		case 14: digraphCopy(graph,g14).nodeRef(nm14).run();
		case 13: digraphCopy(graph,g13).nodeRef(nm13).run();
		case 12: digraphCopy(graph,g12).nodeRef(nm12).run();
		case 11: digraphCopy(graph,g11).nodeRef(nm11).run();
		case 10: digraphCopy(graph,g10).nodeRef(nm10).run();
		case 9: digraphCopy(graph,g9).nodeRef(nm9).run();
		case 8: digraphCopy(graph,g8).nodeRef(nm8).run();
		case 7: digraphCopy(graph,g7).nodeRef(nm7).run();
		case 6: digraphCopy(graph,g6).nodeRef(nm6).run();
		case 5: digraphCopy(graph,g5).nodeRef(nm5).run();
		case 4: digraphCopy(graph,g4).nodeRef(nm4).run();
		case 3: digraphCopy(graph,g3).nodeRef(nm3).run();
		case 2: digraphCopy(graph,g2).nodeRef(nm2).run();
		case 1: digraphCopy(graph,g1).nodeRef(nm1).run();
	}	

	cout << "Finished making copies" << endl;

	//for (set_of_pairs::const_iterator itr = safe_pairs.begin(); itr != safe_pairs.end(); ++itr) 
	omp_set_dynamic(0);
	#pragma omp parallel for num_threads(N_THREADS) shared(ret)
	for (int i = 0; i < n_safe_pairs; i++)
	{

		//************** progress ************** //
		#pragma omp critical
		{
			progress = progress + 1;
		}

		if (progress % 1000 == 0)
		{
			cout << "Safe pair #" << progress << " / " << n_safe_pairs << " ";
		  	cout << "Time: " << currentDateTime();
		  	cout.flush();
		}
		//************** progress ************** //

		StaticDigraph::Arc first_arc = graph.arc(safe_pairs_v[i].first);
		StaticDigraph::Arc second_arc = graph.arc(safe_pairs_v[i].second);
		
		list<int> current_walk;
		current_walk.push_back(graph.id(graph.source(first_arc)));
		current_walk.push_back(graph.id(graph.source(second_arc)));
		current_walk.push_back(graph.id(graph.target(second_arc)));

		// cout << graph.id(graph.source(first_arc)) << ", " << graph.id(graph.source(second_arc)) << ", " << graph.id(graph.target(second_arc)) << endl;

		#pragma omp critical
		{
			reported_arcs.insert(safe_pairs_v[i].first);
			reported_arcs.insert(safe_pairs_v[i].second);	
		}
		
		// problematic case of a 2-cycle a,b,a
		// if so, insert a,b,a as contig and continue
		if (graph.source(first_arc) == graph.target(second_arc))
		{
			contig entry;
			entry.nodes = current_walk;
			#pragma omp critical
			{
				ret.push_back(entry);
			}
			continue;
		}


		vector<int> bad_in_nbrs;
		for (StaticDigraph::InArcIt in_arc2(graph,graph.target(first_arc)); in_arc2 != INVALID; ++in_arc2)
		{
			if (in_arc2 != first_arc)
			{
				bad_in_nbrs.push_back(graph.id(graph.source(in_arc2)));
			}
 		}

		int tid = omp_get_thread_num();
		tid++;
    	extend_pair(graph,GET_GRAPH(tid),GET_NODE_MAP(tid), safe_pairs, second_arc, current_walk, bad_in_nbrs, ret, reported_arcs);
	}

	// including also the arcs not reported in some ec-contig
	for (StaticDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (reported_arcs.count(graph.id(a)) == 0)	
		{
			contig entry;
			entry.nodes.push_back(graph.id(graph.source(a)));
			entry.nodes.push_back(graph.id(graph.target(a)));
			ret.push_back(entry);				
		}
	}

	return ret;
}


int main(int argc, char **argv)
{
	StaticDigraph graph;
	StaticDigraph::NodeMap<size_t> length(graph);
	StaticDigraph::NodeMap<size_t> seqStart(graph);
	set_of_pairs safe_pairs(1 << 17,hash_pair);

	string sequence;
	size_t seqLength;
	size_t kmersize;
	string inputFileName;
    string genome_type = "linear";
    bool circular_genome = false;
	bool build_only = false;
	bool input_from_reads = false;
	int abundance = 1;
	N_THREADS = 1;


	// command line argument parser
	string usage = "\n  %prog OPTIONS"
		"\n\nBrief example:"
		"\n  %prog -i <input file without extension> [--fastq] -k 31 -a 2 [-t 4]";
	const string version = "%prog 0.2\nCopyright (C) 2014-2015 Alexandru Tomescu & Paul Medvedev\n"
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

	parser.add_option("-i", "--input") .type("string") .dest("i") .set_default("") .help("input file, without extension (default assumed extension .fa)");
	parser.add_option("-q", "--fastq") .action("store_true") .dest("fastq_input") .help("add this flag if the input is a fastq file");
	parser.add_option("-k", "--kmersize") .type("int") .dest("k") .action("store") .set_default(31) .help("kmer size (default: %default)");
	parser.add_option("-a", "--abundance") .type("int") .dest("a") .action("store") .set_default(1) .help("minimum abundance (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(1) .help("number of threads, in [1..16] (default: %default)");
	parser.add_option("-g", "--genome-type") .action("store") .dest("g") .type("string") .set_default("linear") .help("genome type: linear|circular (default: %default)");
	parser.add_option("-b", "--build-only") .action("store_true") .dest("build_only") .help("build the de bruijn graph and then exit");
	
	optparse::Values& options = parser.parse_args(argc, argv);

	inputFileName = (string) options.get("i");
	input_from_reads = (options.get("fastq_input") ? true : false);
	kmersize = (size_t) options.get("k");
	abundance = (int) options.get("a");
	N_THREADS = (int) options.get("t");
	build_only = (options.get("build_only") ? true : false);
	genome_type = (string) options.get("g");

	if (inputFileName == "")
	{
		cerr << "Parameter -i | --input is required" << endl;
		cerr << "Use --help to see OPTIONS" << endl;
		return EXIT_FAILURE;
	}

	
	if (genome_type.compare("circular") == 0) 
	{
		circular_genome = true;
	}
	else if (genome_type.compare("linear") == 0) 
	{
		circular_genome = false;
	}
	else 
	{
		cerr << "Genome type (option -g) must be either \'linear\' or \'circular\'" << endl;
		return EXIT_FAILURE;
	}

	if ((N_THREADS < 0) or (N_THREADS > 16))
	{
		cerr << "The number of threads must be between 1 and 16" << endl;
		return EXIT_FAILURE;
	}

	if (input_from_reads)
	{
		if (EXIT_SUCCESS != load_data_from_reads(sequence, seqLength, inputFileName, kmersize, circular_genome, abundance, graph, length, seqStart, safe_pairs))
		{
			return EXIT_FAILURE;
		}
	}
	else
	{
		if (EXIT_SUCCESS != load_data(sequence, seqLength, inputFileName, kmersize, circular_genome, graph, length, seqStart, safe_pairs))
		{
			return EXIT_FAILURE;
		}		
	}


	std::cout.setf(std::ios::boolalpha);
	cout << "Graph has " << countNodes(graph) << " nodes" << endl;
	cout << "Graph has " << countArcs(graph) << " arcs" << endl;
	cout << "Graph has " << countStronglyConnectedComponents(graph) << " strongly connected components" << endl;
	cout << "Graph has unary nodes: " << count_unary_nodes(graph) << endl;
	if (parallelFree(graph) == false)
	{
		cerr << "The graph still has parallel arcs. Contact a developer." << endl;
		return EXIT_FAILURE;
	}
	if (count_unary_arcs(graph) != 0)
	{
		cerr << "Graph still has unary arcs. Contact a developer." << endl;	
		return EXIT_FAILURE;
	}

	if (build_only)
	{
		return EXIT_SUCCESS;
	}

	cout << "Time: " << currentDateTime();
	vector<contig> unitigs;
	unitigs = compute_unitigs(graph, length, seqStart, kmersize, sequence, inputFileName);
	populate_with_strings(sequence, kmersize, graph, seqStart, length, unitigs);
	cout << "----------------------" << endl;
	cout << "Statistics on unitigs:" << endl;
	cout << "----------------------" << endl;
	compute_statistics(unitigs, seqLength);
	print_collection(unitigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".unitig");
	cout << "Time: " << currentDateTime();

	vector<contig> nice_contigs;
	nice_contigs = compute_nice_contigs(graph, length, seqStart, kmersize, sequence, inputFileName);
	populate_with_strings(sequence, kmersize, graph, seqStart, length, nice_contigs);
	cout << "----------------------" << endl;
	cout << "Statistics on nice contigs:" << endl;
	cout << "----------------------" << endl;
	compute_statistics(nice_contigs, seqLength);
	print_collection(nice_contigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".nicecontig");
	cout << "Time: " << currentDateTime();

	// vector<contig> ec_contigs_2;
	// ec_contigs_2 = compute_ec_contigs_2(graph, length, seqStart, safe_pairs, kmersize, sequence, inputFileName);
	// populate_with_strings(sequence, kmersize, graph, seqStart, length, ec_contigs_2);
	// cout << "----------------------" << endl;
	// cout << "Statistics on ec_contigs_2:" << endl;
	// cout << "----------------------" << endl;
	// compute_statistics(ec_contigs_2, seqLength);
	// print_collection(ec_contigs_2, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".ec_contigs_2");

	cout << "safe_pairs.size(): " << safe_pairs.size() << endl;
	if (safe_pairs.size() == 0)
	{
		cout << "Will compute safe pairs" << endl;
		vector<contig> ec_contigs_2;
		ec_contigs_2 = compute_ec_contigs_2(graph, length, seqStart, safe_pairs, kmersize, sequence, inputFileName);
		save_safe_pairs(safe_pairs, inputFileName, kmersize, genome_type);
		//populate_with_strings(sequence, kmersize, graph, seqStart, length, ec_contigs_2);
		//print_collection(ec_contigs_2, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".ec_contig_2");
	}

	vector<contig> ec_contigs;
	try_loading_ec_contigs(ec_contigs,inputFileName,kmersize,genome_type);
	if (ec_contigs.size() == 0)
	{
		ec_contigs = compute_ec_contigs(graph, length, seqStart, safe_pairs);	
	}
	
	// remove_non_maximal_contigs(ec_contigs); 
	populate_with_strings(sequence, kmersize, graph, seqStart, length, ec_contigs);
	cout << "----------------------" << endl;
	cout << "Statistics on ec_contigs:" << endl;
	cout << "----------------------" << endl;
	compute_statistics(ec_contigs, seqLength);
	print_collection(ec_contigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".ec_contig");
	cout << "Time: " << currentDateTime();

	return EXIT_SUCCESS;
}