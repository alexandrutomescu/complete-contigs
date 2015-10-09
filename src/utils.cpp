#include "utils.h"

struct tm * timeinfo;

size_t hash_pair( const pair_of_ints& p )
{
    return p.first ^ p.second;
}

typedef unordered_set<pair_of_ints,function<decltype(hash_pair)>> set_of_pairs;


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

// Check if a file is readable
bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 

size_t count_unary_nodes(StaticDigraph& graph)
{
	size_t n_unary_nodes = 0;
	for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) == 1))
		{
			n_unary_nodes++;
		}
	}
	
	return n_unary_nodes;
}

size_t count_unary_arcs(StaticDigraph& graph)
{
	size_t n_unary_arcs = 0;
	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		if ((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.target(arc)) == 1))
		{
			n_unary_arcs++;
		}
	}
	
	return n_unary_arcs;
}

size_t count_unary_arcs_ld(ListDigraph& graph)
{
	size_t n_unary_arcs = 0;
	for (ListDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		if ((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.target(arc)) == 1))
		{
			n_unary_arcs++;
		}
	}
	
	return n_unary_arcs;
}

void populate_with_strings(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const StaticDigraph::NodeMap<size_t>& length,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			StaticDigraph::Node node = graph.node(*itr);
			string s = sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void populate_with_strings_from_node_labels(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			StaticDigraph::Node node = graph.node(*itr);
			string s = nodeLabel[node]; //sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void populate_with_strings_list_digraph(const string& sequence, 
	size_t kmersize, 
	const ListDigraph& graph,
	const ListDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			ListDigraph::Node node = graph.nodeFromId(*itr);
			string s = nodeLabel[node]; //sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void print_collection(vector<contig>& collection, string inputFileName, string suffix) 
{

	ofstream unitigsFile;
	unitigsFile.open(inputFileName + suffix);

	for (size_t i = 0; i < collection.size(); i++)
	{
		unitigsFile << "> ";
		std::list<int>::const_iterator iterator = collection[i].nodes.begin();
		unitigsFile << *iterator;

		for (++iterator; iterator != collection[i].nodes.end(); ++iterator) 
		{
    		unitigsFile << "," << *iterator;
		}
		unitigsFile << endl;
		unitigsFile << collection[i].str << endl;
	}

	unitigsFile.close();

}

void compute_statistics(vector<contig>& collection, size_t seqLength)  
{
	size_t total_size = 0;
	size_t n_strings = collection.size();
	size_t sum_of_squares = 0;
	size_t max_length = 0;
	size_t n_nodes = 0;
	size_t max_n_nodes = 0;

	for (size_t i = 0; i < collection.size(); i++)
	{
		total_size += collection[i].str.length();
		sum_of_squares += collection[i].str.length() * collection[i].str.length();
		if (collection[i].str.length() > max_length)
		{
			max_length = collection[i].str.length();
		}
		n_nodes += collection[i].nodes.size();
		if (collection[i].nodes.size() > max_n_nodes)
		{
			max_n_nodes = collection[i].nodes.size();
		}
	}

	cout << "Number of strings: " << n_strings << endl;
	cout << "Total length of the strings: " << total_size << endl;
	cout << "Average length: " << total_size / (double)n_strings << endl;
	cout << "Maximum length: " << max_length << endl;
	cout << "Average n. nodes: " << n_nodes / (double)n_strings << endl;
	cout << "Maximum n. nodes: " << max_n_nodes << endl;
	cout << "E-size normalized to assembly length: " << sum_of_squares / (double)total_size << endl;
	//cout << "E-size normalized to genome length: " << sum_of_squares / (double)seqLength << endl;

}

void save_safe_pairs(set_of_pairs& safe_pairs, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	ofstream spFile;
	spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
	for (set_of_pairs::iterator itr = safe_pairs.begin(); itr != safe_pairs.end(); ++itr) 
	{
		spFile << itr->first << " " << itr->second << endl;
	}
	spFile.close();
}

void try_loading_omnitigs(vector<contig>& omnitigs,
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	string ecFileName = inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".omnitigs.nonmaximal";
	ifstream ecFile;
}


void print_graph_in_dot(StaticDigraph& graph, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	ofstream dotFile;
	dotFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".dot");
	dotFile << "digraph DB {" << endl;

	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		dotFile << graph.id(graph.source(arc)) << " -> " << graph.id(graph.target(arc)) << ";" << endl;
	}
	
	dotFile << "}" << endl;
	dotFile.close();
}

void contract_arcs(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize)
{
	vector<ListDigraph::Arc> arcsToContract;


	size_t progress = 0;
	size_t n_arcs = countArcs(graph);
	arcsToContract.reserve(n_arcs);

	for (ListDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (progress % 100000 == 0)
		{
			cout << "Compacting phase, arc : #" << progress << "/" << n_arcs << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		ListDigraph::Node u,v;
		u = graph.source(a);
		v = graph.target(a);

		// OLD: if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,v) == 1))
		// NEW: if the arc (u,v) is the 'middle' arc of a unitig, then contract it
		if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,u) == 1) and
			(countOutArcs(graph,v) == 1) and (countInArcs(graph,v) == 1))
		{
			arcsToContract.push_back(a);
		}
	}

	cout << "arcsToContract has size: " << arcsToContract.size() << endl;
	cout << "count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

	for (auto arc : arcsToContract)
	{
		if (progress % 100000 == 0)
		{
			cout << "Compacting phase 2, arc : #" << progress << "/" << arcsToContract.size() << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		if (not graph.valid(arc))
		{
			cout << "arc not valid" << endl;
		}

		ListDigraph::Node u,v;
		u = graph.source(arc);
		v = graph.target(arc);
		length[u] = length[u] + length[v] - (kmersize - 1);

		//graph.contract(u, v);
		for (ListDigraph::OutArcIt out_arc(graph,v); out_arc != INVALID;)
		{
			ListDigraph::OutArcIt next_out_arc = out_arc;
			++next_out_arc;
			graph.changeSource(out_arc,u);
			out_arc = next_out_arc;
		}
		graph.erase(v);

	}

	cout << "after compacting, count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

	cout << "Compacted the graph" << endl;

}

void contract_arcs_from_reads(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	string& sequence)
{
	vector<ListDigraph::Arc> arcsToContract;
	ListDigraph::NodeMap<string> nodeLabels(graph);

	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		nodeLabels[node] = sequence.substr(seqStart[node],length[node]);
	}
	

	size_t progress = 0;
	size_t n_arcs = countArcs(graph);
	arcsToContract.reserve(n_arcs);

	for (ListDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (progress % 100000 == 0)
		{
			cout << "Compacting phase, arc : #" << progress << "/" << n_arcs << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		ListDigraph::Node u,v;
		u = graph.source(a);
		v = graph.target(a);

		// OLD: if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,v) == 1))
		// NEW: if the arc (u,v) is the 'middle' arc of a unitig, then contract it
		if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,u) == 1) and
			(countOutArcs(graph,v) == 1) and (countInArcs(graph,v) == 1))
		{
			arcsToContract.push_back(a);
		}
	}

	cout << "arcsToContract has size: " << arcsToContract.size() << endl;
	cout << "count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

	for (auto arc : arcsToContract)
	{
		if (progress % 100000 == 0)
		{
			cout << "Compacting phase 2, arc : #" << progress << "/" << arcsToContract.size() << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		if (not graph.valid(arc))
		{
			cout << "arc not valid" << endl;
		}

		ListDigraph::Node u,v;
		u = graph.source(arc);
		v = graph.target(arc);
		//length[u] = length[u] + length[v] - (kmersize - 1);
		nodeLabels[u] = plus_strings(nodeLabels[u], nodeLabels[v], kmersize);

		// the following code is doing "graph.contract(u, v);" (this lemon function doesn't behave as expected)
		for (ListDigraph::OutArcIt out_arc(graph,v); out_arc != INVALID;)
		{
			ListDigraph::OutArcIt next_out_arc = out_arc;
			++next_out_arc;
			graph.changeSource(out_arc,u);
			out_arc = next_out_arc;
		}
		graph.erase(v);

	}

	cout << "after compacting, count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

	// creating the new sequence and changing
	// length[] and seqStart[] maps to be relative to this new sequence
	sequence = "";
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		seqStart[node] = sequence.length();
		length[node] = nodeLabels[node].length();
		sequence += nodeLabels[node];
	}

	cout << "Compacted the graph" << endl;

}


void construct_graph(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize, 
	const bool circular_genome,
	const string& sequence,
	const size_t seqLength,
	const int abundance)
{
	unordered_map<std::string,int> kmersToNodes;
	unordered_map<int,int> nodeAbundance;
	unordered_map<std::string,unordered_set<int>> in_nbrs;

	ListDigraph::Node currentNode, previousNode = INVALID;
	std::string currentKmer;
	// std::string currentKmer = sequence.substr(0,kmersize);
	// // cout << "current kmer is: " << currentKmer << endl;
	// previousNode = graph.addNode();
	// kmersToNodes[currentKmer] = graph.id(previousNode);
	// length[previousNode] = kmersize;
	// seqStart[previousNode] = 0;

	size_t kmerLimit;
	if (circular_genome)
	{
		kmerLimit = seqLength;
	} else
	{
		kmerLimit = seqLength - kmersize + 1;
	}

	size_t progress = 0;

	// traversing the sequence
	for (size_t i = 0; i <= kmerLimit; i++)
	{

		if (progress % 1000000 == 0)
		{
			cout << "Sequence position: #" << progress << "/" << kmerLimit << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		// if new kmer, then add a node to the graph
		currentKmer = sequence.substr(i,kmersize);

		// if current kmer contains # then skip it
		if (currentKmer.find("#") != std::string::npos)
		{
			previousNode = INVALID;
		}
		else 
		{
			// current kmer does not contain #
			// cout << "current kmer is: " << currentKmer << endl;
			if (kmersToNodes.count(currentKmer) > 0) 
			{
				currentNode = graph.nodeFromId(kmersToNodes[currentKmer]);
				nodeAbundance[kmersToNodes[currentKmer]]++;
			}
			else
			{
				currentNode = graph.addNode();
				//kmersToNodes.insert(std::make_pair<std::string,StaticDigraph::Node>(currentKmer, currentNode));
				kmersToNodes[currentKmer] = graph.id(currentNode);
				nodeAbundance[kmersToNodes[currentKmer]] = 1;
				length[currentNode] = kmersize;
				seqStart[currentNode] = i;
				// cout << "inserted seq start" << seqStart[kmersToNodes[currentKmer]] << endl;
			}

			if (previousNode != INVALID)
			{
				if (in_nbrs[currentKmer].count(graph.id(previousNode)) == 0)
				{
					graph.addArc(previousNode, currentNode);
					in_nbrs[currentKmer].insert(graph.id(previousNode));
				}	
			}
			previousNode = currentNode;
		}
	}	

	// traversing the graph and removing the nodes with abundance strictly lower than abundance
	vector<ListDigraph::Node> nodesToRemove;

	// finding the nodes to remove
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if (nodeAbundance[graph.id(node)] < abundance)
		{
			nodesToRemove.push_back(node);
		}
	}

	// removing the nodes from nodesToRemove
	for (auto node : nodesToRemove)
	{
		graph.erase(node);
	}


	cout << "Constructed the graph" << endl;

}


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
)
{
	string genome_type;
	if (circular_genome)
	{
		genome_type = "circular";
	} else
	{
		genome_type = "linear";
	}

	try 
	{
		ifstream sequenceFile;
		string line;
        sequenceFile.open(inputFileName);
        getline(sequenceFile, line); // reading the header
        while (getline(sequenceFile, line))
        {
        	sequence += line;
        }

        seqLength = sequence.length();
        make_upper_case(sequence);

        if (circular_genome)
        {
        	sequence = sequence + sequence; // + sequence.substr(0,kmersize);
        }
        sequenceFile.close();  
	} catch (Exception& error) 
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

	if (! is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf"))
	{
		ListDigraph temporary_graph;
		ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
		ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

		construct_graph(temporary_graph, temporary_length, temporary_seqStart, kmersize, circular_genome, sequence, seqLength, 1);
		if (not do_not_contract_arcs)
		{
			contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		}
		//contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		// saving graph to file
		try 
		{
			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			digraphWriter(graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", length).
				nodeMap("seqStart", seqStart).
				run();
			// print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	else
	{
		try 
		{
			ListDigraph temporary_graph;
			ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
			ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

			digraphReader(temporary_graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", temporary_length).
				nodeMap("seqStart", temporary_seqStart).
				run();

			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			// attach the node labels
			for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
			{
				nodeLabel[node] = sequence.substr(seqStart[node],length[node]);
			}


		}
		catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs"))
	{
		int first, second;
		ifstream spFile;
		spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
		while (spFile >> first)
		{
			spFile >> second;
			safe_pairs.insert(pair_of_ints(first,second));
		}
		spFile.close();
		cout << safe_pairs.size() << endl;
	}

	return EXIT_SUCCESS;
}

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
)
{
	string genome_type;
	if (circular_genome)
	{
		genome_type = "circular";
	} else
	{
		genome_type = "linear";
	}

	if (! is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf"))
	{
		try 
		{
			ifstream sequenceFile;
			string line;
	        sequenceFile.open(inputFileName);
	        
	        while (getline(sequenceFile, line)) // reading the header
	        {
	        	getline(sequenceFile, line); // the read itself
	        	make_upper_case(line);
	        	sequence += line + "#";
	        	sequence += reverse_complement(line) + "#";
	        	getline(sequenceFile, line);
	        	getline(sequenceFile, line);
	        }

	        seqLength = sequence.length();
	        make_upper_case(sequence);

	        sequenceFile.close();  
		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}

		ListDigraph temporary_graph;
		ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
		ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

		construct_graph(temporary_graph, temporary_length, temporary_seqStart, kmersize, circular_genome, sequence, seqLength, abundance);
		if (not do_not_contract_arcs)
		{
			contract_arcs_from_reads(temporary_graph, temporary_length, temporary_seqStart, kmersize, sequence);	
		}
		//contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		// saving graph to file
		try 
		{
			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			digraphWriter(graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", length).
				nodeMap("seqStart", seqStart).
				attribute("sequence", sequence).
				run();
			// print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	else
	{
		try 
		{
			ListDigraph temporary_graph;
			ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
			ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

			digraphReader(temporary_graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", temporary_length).
				nodeMap("seqStart", temporary_seqStart).
				attribute("sequence", sequence).
				run();

			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}

			print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		}
		catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs"))
	{
		int first, second;
		ifstream spFile;
		spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
		while (spFile >> first)
		{
			spFile >> second;
			safe_pairs.insert(pair_of_ints(first,second));
		}
		spFile.close();
		cout << safe_pairs.size() << endl;
	}

	return EXIT_SUCCESS;
}

