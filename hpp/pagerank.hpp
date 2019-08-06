#ifndef _PAGERANK_H_
#define _PAGERANK_H_

#include <iostream>
#include <vector>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <chrono>

using namespace std;

class pagerank{
private:
	const int gauss_start = 1;
	double alpha = 0.85;
	double epsi;
	int xNum;
	int t;
	int clusterNum;

	void normalize_outedge_weight(graph& g);
	void pagerank_from_scratch(graph& g);
	void get_rank_for_rankclus(graph& g);

public:
	pagerank(int _xNum, int _t, int clusterNum, double _epsi);
	void solve(graph& g, int iteration_num);
};

#endif 
