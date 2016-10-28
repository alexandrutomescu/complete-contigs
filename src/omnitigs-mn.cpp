#include "utils.h"

// O(nm) time computation of omnitigs
// TODO: Omnitigs containing non-strong-bridges are not computed, yet

class OmnitigNM {
	StaticDigraph& G;

	OutDegMap<StaticDigraph> out_deg;
	InDegMap<StaticDigraph> in_deg;

	StaticDigraph::ArcMap<bool> redundant; // not a strong bridge

	/*
	 * prec[e] == f if $f \prec e$,
	 * and prec[e] == INVALID otherwise
	 */
	StaticDigraph::ArcMap<StaticDigraph::Arc> prec;

	StaticDigraph::ArcMap<bool> compiling_started;

	StaticDigraph::ArcMap<Path<StaticDigraph>> omnitig_ending_with;

	/*
	 * longest_suffix_length[e] = length of the longest suffix of pe
	 * such that there are no back-paths from a sibling branch e'' of e'
	 */
	StaticDigraph::ArcMap<int> longest_suffix_length;

	/*
	 * whether the left-maximal omnitig $wep$,
	 * where $p$ is univocal,
	 * is also right-maximal.
	 */
	StaticDigraph::ArcMap<bool> maximal;

public:
	OmnitigNM(StaticDigraph& G) :
					G(G),
					out_deg(G),
					in_deg(G),
					redundant(G, true),
					prec(G, INVALID),
					compiling_started(G, false),
					omnitig_ending_with(G),
					longest_suffix_length(G, -1),
					maximal(G, true) {
	}

	vector<contig> run() {
		cout << "Entered omnitigs O(nm)" << endl;

		cout << "Computing branches/strong bridge candidates...";
		compute_strong_bridge_candidates();
		cout << " done." << endl;

		return compute_all_omnitigs();
	}

private:
	StaticDigraph::Node extend_backward_univocal(
			StaticDigraph::Node node,
			Path<StaticDigraph>& path) {
		assert(path.length() == 0 || pathSource(G, path) == node);

		while (in_deg[node] == 1) {
			StaticDigraph::Arc next = StaticDigraph::InArcIt(G, node);
			path.addFront(next);
			node = G.source(next);
		}

		assert(checkPath(G, path));
		assert(path.length() == 0 || pathSource(G, path) == node);

		return node;
	}

StaticDigraph::Node extend_forward_univocal(
			StaticDigraph::Node node,
			Path<StaticDigraph>& path) {
		assert(path.length() == 0 || pathTarget(G, path) == node);

		while (out_deg[node] == 1) {
			StaticDigraph::Arc next = StaticDigraph::OutArcIt(G, node);
			path.addBack(next);
			node = G.target(next);
		}

		assert(checkPath(G, path));
		assert(path.length() == 0 || pathTarget(G, path) == node);

		return node;
	}

void compute_strong_bridge_candidates() {
		StaticDigraph::NodeIt u(G);
		ReverseDigraph<StaticDigraph> GR(G);

		Bfs<StaticDigraph> visit(G);
		visit.run(u);
		Bfs<ReverseDigraph<StaticDigraph>> rvisit(GR);
		rvisit.run(u);

		for (StaticDigraph::NodeIt u(G); u != INVALID; ++u) {
			auto e = visit.predArc(u);
			if (e != INVALID) redundant[e] = false;
			auto f = rvisit.predArc(u);
			if (f != INVALID) redundant[f] = false;
		}
	}

	bool try_closed_path_omnitig(
			StaticDigraph::Arc e,
			const Path<StaticDigraph>& q0) {
		StaticDigraph::Node se = G.source(e);
		if (out_deg[se] > 2) return false;

		StaticDigraph::OutArcIt e1(G, se);
		if (e1 == e) ++e1;

		Path<StaticDigraph> p;
		p.addBack(e1);
		StaticDigraph::Node tp = extend_forward_univocal(G.target(e1), p);

		if (tp != se) return false;

		// p0 is a closed path from s(e)
		assert(checkPath(G, p));
		assert(pathSource(G, p) == se);
		assert(pathTarget(G, p) == se);

		Path<StaticDigraph> w(q0);
		for (Path<StaticDigraph>::ArcIt arc(p); arc != INVALID; ++arc) {
			w.addBack(arc);
		}
		w.addBack(e);

		assert(checkPath(G, w));
		assert(w.length() == q0.length() + p.length() + 1);

		omnitig_ending_with[e] = w;

		return true;
	}
	bool try_redundant_arc_omnitig(
			StaticDigraph::Arc e,
			const Path<StaticDigraph>& q0) {
		if (!redundant[e]) return false;

		Path<StaticDigraph> w = q0;
		w.addBack(e);
		assert(checkPath(G, w));
		omnitig_ending_with[e] = w;

		return true;
	}
Path<StaticDigraph> longest_suffix(
			const Path<StaticDigraph>& w,
			StaticDigraph::Arc e,
			const Bfs<FilterArcs<StaticDigraph> >& visit) {
		Path<StaticDigraph> suffix;
		suffix.addBack(e);
		for (int i = w.length() - 1; i >= 0; i--) {
			StaticDigraph::Arc f = w.nth(i);
			bool stop = false;
			for (StaticDigraph::InArcIt f1(G, G.target(f)); f1 != INVALID;
					++f1) {
				if (f1 == f || f1 == e) continue;

				if (visit.reached(G.source(f1))) {
					stop = true;
					break;
				}
			}
			if (stop) break;

			suffix.addFront(f);
		}
		assert(checkPath(G, suffix));
		return suffix;
	}
void process_branch(StaticDigraph::Arc e) {
		StaticDigraph::Node se = G.source(e);

		assert(out_deg[se] >= 2);

		Path<StaticDigraph> q0;
		StaticDigraph::Node sq0 = se;
		sq0 = extend_backward_univocal(sq0, q0);

		if (try_closed_path_omnitig(e, q0)) return;
		if (try_redundant_arc_omnitig(e, q0)) return;

		StaticDigraph::ArcMap<bool> Ge_map(G, true);
		FilterArcs<StaticDigraph> Ge(G, Ge_map);
		Ge.disable(e);

		Bfs<FilterArcs<StaticDigraph>> visit(Ge);
		visit.init();
		visit.addSource(se);

		// Optimization: visit at most two in-neighbors of s(q0)
		FilterArcs<StaticDigraph>::NodeMap<bool> target(Ge, false);
		for (StaticDigraph::InArcIt g(G, sq0); g != INVALID; ++g) {
			target[G.source(g)] = true;
		}
		visit.start(target);
		visit.start(target);

		// make a closed path
		StaticDigraph::InArcIt g(G, sq0);
		for (; g != INVALID; ++g) {
			if (visit.reached(G.source(g))) break;
		}
		assert(g != INVALID);

		Path<StaticDigraph> e1p(visit.path(G.source(g)));
		e1p.addBack(g);
		for (int i = 0; i < q0.length(); i++) {
			e1p.addBack(q0.nth(i));
		}

		// e1pq0 is a closed path from s(e)
		assert(checkPath(G, e1p));
		assert(pathSource(G, e1p) == G.source(e));
		assert(pathTarget(G, e1p) == G.source(e));

		Path<StaticDigraph> suffix = longest_suffix(e1p, e, visit);
		for (int i = suffix.length() - 2; i >= 0; i--) {
			StaticDigraph::Arc f = suffix.nth(i);

			if (out_deg[G.source(f)] >= 2) {
				assert(prec[e] == INVALID);
				assert(f != e);
				prec[e] = f;
				break;
			}
		}

		if (prec[e] == INVALID) {
			omnitig_ending_with[e] = suffix;
		} else {
			longest_suffix_length[e] = suffix.length();
		}
	}

	void compile_omnitig_ending_with(StaticDigraph::Arc e) {
		StaticDigraph::Arc f = prec[e];
		if (f == INVALID) return;
		if (omnitig_ending_with[e].length()) return;

		assert(!compiling_started[e]);
		compiling_started[e] = true;

		compile_omnitig_ending_with(f);

		Path<StaticDigraph> w = omnitig_ending_with[f];
		extend_forward_univocal(pathTarget(G, w), w);
		w.addBack(e);

		assert(checkPath(G, w));

		assert(w.length() >= longest_suffix_length[e]);
		while (w.length() > longest_suffix_length[e]) {
			w.eraseFront();
		}

		omnitig_ending_with[e] = w;
	}

	vector<contig> compute_all_omnitigs() {
		for (StaticDigraph::ArcIt e(G); e != INVALID; ++e) {
			if (out_deg[G.source(e)] == 1) continue;
			process_branch(e);
		}

		vector<contig> ret;

		for (StaticDigraph::ArcIt e(G); e != INVALID; ++e) {
			if (out_deg[G.source(e)] == 1) continue;
			compile_omnitig_ending_with(e);

			Path<StaticDigraph> w = omnitig_ending_with[e];
			extend_forward_univocal(G.target(e), w);

			assert(checkPath(G, w));

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

vector<contig> compute_omnitigs_nm(StaticDigraph& G) {
	return OmnitigNM(G).run();
}

