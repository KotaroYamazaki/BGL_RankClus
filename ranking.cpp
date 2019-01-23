#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

// 行列ベクトル積
graph ranking(graph &g){
    vertex_iterator i, j;
    //iterator を用いて最初のiterから最後まで回す
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        double tmp = 0;
        //出る全てのエッジに対してループを回す
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            // ノード *i の入エッジの重み（g[*e.first].weight）と
            // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
            tmp += g[*e.first].weight * g[source(*e.first, g)].previous_rank;
        }
        g[*i].next_rank = tmp;
    }
    return g;
}