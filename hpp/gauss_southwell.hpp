#ifndef _GAUSS_SOUTHWELL_H_
#define _GAUSS_SOUTHWELL_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <queue>

using namespace std;

class gauss_southwell
{
private:
	double epsi;
	double alpha;
	int clusterNum;

	void calc_initial_residual(graph &g);
	void init_pregraph(graph& g);
	pair<queue<int>, vector<bool>> calc_tracking_residual(graph &g, int clusterNum);


public:
	gauss_southwell(double _epsi, double _alpha, int _clusterNum);
	void init(graph &g);
	void update_pregraph(graph &g, int index);
	void solve(graph &g);
};

#endif