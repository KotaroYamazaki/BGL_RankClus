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
extern int xNum;

graph clustering(graph &g, vector<graph>& subgraph){
    vector<vector<double>> pi;
    vector<double> p;
    vertex_iterator i,j;
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        p.push_back(1.0*WkXY_sum[clusterNum]/WXY_sum);
    }

    for(int z = 0; z < K; z++){
        
        double tmp_p = 0;
        for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                //サブグラフでエッジないとこにあくせすしてしまう？　
                //それは違うか
                tmp_p += g[*e.first].weight * subgraph[z][target(*e.first, g)].rx *subgraph[z][source(*e.first, g)].ry * p[z];
            }
        }
        p[z] = tmp_p/WXY_sum;
    }

    // calc pi
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
        double tmp_sum = 0;
        for(int l = 0; l < K; l++){
            tmp_sum = subgraph[l][*i].rx * p[l];
        }

        for(int z = 0; z < K; z++){ 
            pi[z][*i] = subgraph[z][*i].conditional_rank * p[z];
            pi[z][*i] /= tmp_sum;
        }
    }

    vector<vector<double>> s(xNum);
    for(int i = 0; i < xNum; i++){
        for(int k = 0; k < K; k++){
            s[i].push_back(pi[i][k]);
        }
    }
    
    
    return g;
}