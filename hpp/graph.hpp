#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>

using namespace std;

enum vertex_trajectory{
    NOTHING = 0,
    STAY = 1,
    ADD = 2,
    LEAVE = 3,
};

// ノードの属性値を自由に定義
class vertex_property
{
    public:
        string label;
        string name;
        int int_descriptor;
        //ランク情報
        double ry;
        double rx;
        double p_rank;
        double conditional_rank;
        vector<vertex_trajectory> state;
        // クラスタ所属情報
        int cluster_label;
        bool same_previous_cluster;

        bool belongs_to_cluster(int _label){
            return (cluster_label == _label);
        }
};

// エッジの属性値を自由に定義
class edge_property
{
    public:
        string label;
        double weight;
};

// グラフの属性値を自由に定義
class graph_property
{
    //アクセス方法がわからない
    public:
        int _xNum;
        int _yNum;
        double _edge_sum;
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
    //boost::undirectedS,
    boost::bidirectionalS,
    vertex_property,
    edge_property,
    graph_property
    >;

// 簡単のために、名前をつける
using edge = pair<int, int>;
using vertex_iterator = boost::graph_traits<graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<graph>::edge_iterator;
using vertex_descriptor = boost::graph_traits<graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<graph>::edge_descriptor;
using adjacency_iterator = boost::graph_traits<graph>::adjacency_iterator;

#endif
