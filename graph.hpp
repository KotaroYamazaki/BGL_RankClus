#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>

// ノードの属性値を自由に定義
struct vertex_property
{
    std::string label;
    double previous_rank;
    double next_rank;
    int int_descriptor;
};

// エッジの属性値を自由に定義
struct edge_property
{
    std::string label;
    double weight;
};

// グラフの属性値を自由に定義
struct graph_property
{
};

// グラフ構造を定義
// CSR形式で保持 -> compressed_sparse_row_graph
// 隣接リスト -> adjacency_list
//
// directed graph -> boost::directedS
// undirected graph -> boost::undirectedS
// bidirected graph -> boost::bidirectionalS
// あるノードから in_edge と out_edge の両方たどりたい場合は、bidirectionalS が良い
using graph = boost::compressed_sparse_row_graph<
    boost::bidirectionalS,
    vertex_property,
    edge_property,
    graph_property
    >;

// 簡単のために、名前をつける
using edge = std::pair<int, int>;
using vertex_iterator = boost::graph_traits<graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<graph>::edge_iterator;
using vertex_descriptor = boost::graph_traits<graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<graph>::edge_descriptor;
using adjacency_iterator = boost::graph_traits<graph>::adjacency_iterator;

#endif
