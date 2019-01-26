#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;
extern vector<int> WkXY_sum;
extern int K;
extern int WXY_sum;

graph clustering(graph &g, vector<graph>& subgraph){
    vector<vector<double>> pi;
    vector<double> p;
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        p.push_back(1.0*WkXY_sum[clusterNum]/WXY_sum);
    }
    return g;
}