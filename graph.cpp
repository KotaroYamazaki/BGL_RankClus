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

int main()
{
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    // グラフの詳細を出力
	print_detail(g);

    // 行列ベクトル積
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        double tmp = 0;
        
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            // ノード *i の入エッジの重み（g[*e.first].weight）と
            // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
            tmp += g[*e.first].weight * g[source(*e.first, g)].previous_rank;
        }

        g[*i].next_rank = tmp;
    }
}

struct edge_property add_edge_property(string s, int w){
    struct edge_property a;
    a.label = s;
    a.weight = w;
    return a;
}

graph construct_graph(){
    // エッジのリスト
    std::vector<edge> edge_vector;
    // 各エッジの属性値の構造体のリスト
    std::vector<edge_property> property_vector;

    string str;
    string fnameWXY = "AVWeight.csv";
    string fnameWYY = "AAWeight.csv";
	ifstream ifs(fnameWXY);
    int from, to, val;
    // Target type
	while(getline(ifs,str)){
		sscanf(str.c_str(), "%d %d %d", &from, &to, &val);
        struct edge_property a;
        a.label = "target";
        a.weight = val;
        property_vector.push_back(a);
		edge_vector.push_back(edge(from, to));
	}
    //Attribute Type
    ifstream ifs2(fnameWYY);
    while(getline(ifs2,str)){
		sscanf(str.c_str(), "%d %d %d", &from, &to, &val);
        struct edge_property a;
        a.label = "attribute";
        a.weight = val;
        property_vector.push_back(a);
		edge_vector.push_back(edge(from, to));
	}

    // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
    auto tag = boost::edges_are_unsorted_multi_pass;
    // グラフのコンストラクタ
    // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
    graph g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), 5639);

    return g;
}

// graph construct_graph()
// {
//     // エッジのリスト
//     std::vector<edge> edge_vector;
//     // 各エッジの属性値の構造体のリスト
//     std::vector<edge_property> property_vector;


//     int from, to;
//     // 0 -> 1
//     //struct edge_property a;
//     from = 0, to = 1;
//     // a.label = "target";
//     // a.weight = 1.0;
//     struct edge_property a = add_edge_property("target", 1);
//     property_vector.push_back(a);
//     edge_vector.push_back(edge(from, to));

//     // 0 -> 2
//     struct edge_property b;
//     from = 0, to = 2;
//     b.label = "target";
//     b.weight = 1.0;
//     property_vector.push_back(b);
//     edge_vector.push_back(edge(from, to));

//     // 2 -> 1
//     struct edge_property c;
//     from = 2, to = 1;
//     c.label = "attribute";
//     c.weight = 1.0;
//     property_vector.push_back(c);
//     edge_vector.push_back(edge(from, to));


//     // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
//     auto tag = boost::edges_are_unsorted_multi_pass;
//     // グラフのコンストラクタ
//     // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
//     graph g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), 3);

//     return g;
// }

// グラフの初期化
void init_graph(graph& g)
{
    // ノードの走査
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        // (*vertex_iterator) で vertex_descriptor になる
        // g[vertex_descriptor].属性（vertex_property で定義したもの）で参照できる
        g[*i].label = "target";
        //g[*i].previous_rank = 1.0/num_vertices(g);
        g[*i].next_rank = 1.0/num_vertices(g);
        g[*i].int_descriptor = static_cast<int>(*i);
    }
}

void print_detail(graph& g)
{
    cout << "# vertices: " << num_vertices(g) << endl;
    cout << "# edges   : " << num_edges(g) << endl;

    vertex_iterator i, j;
    cout << " --- adjacent vertices --- " << endl;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        cout <<  *i << " --> " << flush;
        // adjacent_vertices(vertex_descriptor, graph) で隣接ノードを取得できる
        for (auto itr = adjacent_vertices(*i, g); itr.first!=itr.second; itr.first++) {
            cout << *itr.first << ", " << flush;
        }
        cout << endl;
    }


    // out_edges(vertex_descriptor, graph) でノードの出エッジが取得できる
    cout << " --- edges --- " << endl;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
            // target(edge_descriptor, graph) でエッジの先のノードを取得できる
            // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
            // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
            cout << *i << " --> " << target(*e.first, g) << ": " << g[*e.first].weight << " (" << g[*e.first].label << ")" << endl;
        }
    }
}


