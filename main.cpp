#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

graph construct_graph();
void init_graph(graph& g);
void print_detail(graph& g);
graph ranking(graph &g);


int main()
{
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    // グラフの詳細を出力
	print_detail(g);

    ranking(g);

}