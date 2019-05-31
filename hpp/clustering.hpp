#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

class clustering{
    private:
    vector<double> WkXY_sum;
    int K;
    double WXY_sum;
    int xNum;
    vector<vector<int>> cluster_label;
    double Norm(vector<double>& array );

    public:
    clustering(vector<double>& _WkXY_sum, int _K, double _WXY_sum,int _xNum, vector<vector<int>>& _cluster_label);
    void begin_clustering(graph &g, vector<graph>& subgraph);
};