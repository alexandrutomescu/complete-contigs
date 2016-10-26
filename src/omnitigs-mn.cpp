#include "utils.h"

// O(nm) time computation of omnitigs
// TODO: Omnitigs containing non-strong-bridges are not computed, yet

class OmnitigNM {
	StaticDigraph& graph;
	StaticDigraph::NodeMap<bool> is_branching;
	StaticDigraph::ArcMap<bool> is_branch;
	StaticDigraph::ArcMap<bool> is_strong_bridge;
	StaticDigraph::ArcMap<bool> computation_started;
	StaticDigraph::ArcMap<Path<StaticDigraph>> omnitig_ending_with;
	vector<int> strong_branches;
	int started_count = 0;

public:
	OmnitigNM(StaticDigraph& G) :
		graph(G),
		is_branching(G, false),
		is_branch(G, false),
		is_strong_bridge(G, false),
		computation_started(G, false),
		omnitig_ending_with(G)
	{}

	vector<contig> run()
	{
		cout << "Entered omnitigs O(nm)" << endl;

		cout << "Computing branches/strong bridge candidates...";
		compute_strong_bridge_candidates();
		compute_branches();
		compute_strong_branches();
		cout << " done." << endl;

		return compute_all_omnitigs();
	}

private:
	void compute_strong_bridge_candidates() {
		auto& G = graph;

		StaticDigraph::NodeIt u(G);
		ReverseDigraph<StaticDigraph> GR(G);

		Dfs<StaticDigraph> visit(G);
		visit.run(u);
		Dfs<ReverseDigraph<StaticDigraph>> rvisit(GR);
		rvisit.run(u);

		for (StaticDigraph::NodeIt u(G); u != INVALID; ++u) {
			auto e = visit.predArc(u);
			if (e != INVALID) is_strong_bridge[e] = true;
			auto f = rvisit.predArc(u);
			if (f != INVALID) is_strong_bridge[f] = true;
		}
	}

	void compute_branches() {
		auto& G = graph;

		for (StaticDigraph::NodeIt u(G); u != INVALID; ++u) {
			int out_deg = countOutArcs(G, u);
			assert(out_deg > 0);
			is_branching[u] = (out_deg > 1);
			if (!is_branching[u]) {
				auto e = StaticDigraph::OutArcIt(G, u);
				is_branch[e] = false;
			} else {
				for (StaticDigraph::OutArcIt e(G, u); e != INVALID; ++e) {
					is_branch[e] = true;
				}
			}
		}
	}

	void compute_strong_branches() {
		auto& G = graph;
		for (StaticDigraph::ArcIt e(G); e != INVALID; ++e) {
			if (is_branch[e] && is_strong_bridge[e]) {
				strong_branches.push_back(G.id(e));
			}
		}
	}

	Path<StaticDigraph> longest_suffix(
			Path<StaticDigraph> const& w,
			StaticDigraph::Arc const& e,
			Dfs<FilterArcs<StaticDigraph>>& visit) {
		auto& G = graph;

		assert(checkPath(G, w));
		assert(pathTarget(G, w) == G.source(e));
		assert(is_strong_bridge[e]);

		int start = 0;
		for(int i = w.length()-1; i >= 0; i--) {
			auto f = w.nth(i);
			for(StaticDigraph::InArcIt f1(G, G.target(f)); f1 != INVALID; ++f1) {
				if(f1 == f || f1 == e) continue;
				if(visit.reached(G.source(f1))) {
					start = i+1;
					break;
				}
			}
		}

		Path<StaticDigraph> ret;

		for(int i = start; i < w.length(); i++) {
			ret.addBack(w.nth(i));
		}
		ret.addBack(e);

		assert(checkPath(G, ret));
		assert(ret.back() == e);

		return ret;
	}

	void compute_omnitig_ending_with(StaticDigraph::Arc const& e)
	{
		auto& G = graph;
		auto& res = omnitig_ending_with[e];

		assert(is_branch[e]);
		assert(is_strong_bridge[e]);

		if (res.length()) {
			assert(computation_started[e]);
			return;
		}
		assert(!computation_started[e]); // check acyclic

		computation_started[e] = true;

		computation_started[e] = true;
		if(++started_count % 1000 == 0) {
			cout << "Strong branch #" << started_count << "/" << strong_branches.size();
			cout << " Time: " << currentDateTime();
		}

		StaticDigraph::ArcMap<bool> Ge_map(G, true);
		FilterArcs<StaticDigraph> Ge(G, Ge_map);

		Ge.disable(e);

		Dfs<FilterArcs<StaticDigraph>> visit(Ge);

		visit.run(G.source(e));

		// find any last edge g which closes a path
		StaticDigraph::InArcIt g(G, G.source(e));
		for(; g != INVALID; ++g) {
			if(visit.reached(G.source(g))) break;
		}
		assert(g != INVALID);

		StaticDigraph::Arc f;
		Path<StaticDigraph> e1p(visit.path(G.source(g)));
		e1p.addBack(g);

		// e1p is a closed path from s(e)
		assert(checkPath(G, e1p));
		assert(pathSource(G, e1p) == G.source(e));
		assert(pathTarget(G, e1p) == G.source(e));

		// first edge of e1p is e1 != e
		assert(e1p.front() != e);

		Path<StaticDigraph> q;
		for(int i = e1p.length()-1; i >= 0; i--) {
			StaticDigraph::Arc fi = e1p.nth(i);
			if (is_branch[fi]) {
				f = fi;
				break;
			}
			q.addFront(fi);
		}

		assert(checkPath(G, q));
		assert(is_branch[f]);

		Path<StaticDigraph> fq(q);
		fq.addFront(f);
		assert(checkPath(G, fq));

		auto w1 = longest_suffix(fq, e, visit);

		assert(checkPath(G, w1));
		assert(w1.back() == e);

		if(w1.length() <= fq.length()) {
			res = w1;
		} else {
			assert(w1.length() == fq.length() + 1);
			assert(w1.front() == f);

			compute_omnitig_ending_with(f);

			Path<StaticDigraph> w2q = omnitig_ending_with[f];
			for(Path<StaticDigraph>::ArcIt g(q); g != INVALID; ++g) {
				w2q.addBack(g);
			}
			assert(checkPath(G, w2q));

			res = longest_suffix(w2q, e, visit);
		}

		assert(checkPath(G, res));
	}

	vector<contig> compute_all_omnitigs()
	{
		auto& G = graph;

		vector<contig> ret;
		for (auto it = strong_branches.cbegin(); it != strong_branches.cend(); ++it) {
			const auto& e = G.arc(*it);
			assert(is_strong_bridge[e]);
			assert(is_branch[e]);

			compute_omnitig_ending_with(e);

			auto w = omnitig_ending_with[e];

			while (!is_branching[pathTarget(G, w)]) {
				w.addBack(StaticDigraph::OutArcIt(G, pathTarget(G, w)));
			}

			contig entry;
			for (Path<StaticDigraph>::ArcIt f(w); f != INVALID; ++f) {
				entry.nodes.push_back(G.id(G.source(f)));
			}
			entry.nodes.push_back(G.id(pathTarget(G, w)));
			ret.push_back(entry);
		}
		return ret;
	}

};

vector<contig> compute_omnitigs_nm(StaticDigraph& G)
{
	return OmnitigNM(G).run();
}

