#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

const int N = 20;
int Wsum = 0;
// struct edge_property add_edge_property(string s, int w){
//     struct edge_property a;
//     a.label = s;
//     a.weight = w;
//     return a;
// }
vector<string> name_vector;

graph construct_graph(){
    // エッジのリスト
    std::vector<edge> edge_vector;
    // 各エッジの属性値の構造体のリスト
    std::vector<edge_property> property_vector;

    string str, name_X, name_Y;
    string fnameWXY = "AVWeight.csv";
    string fnameWYY = "AAWeight2.csv";
    string file_X = "VenueListSimple.txt";
    string file_Y = "AuthorList.txt";
    
	ifstream ifs_WXY(fnameWXY);
    ifstream ifs_X(file_X);
    ifstream ifs_WYY(fnameWYY);
    ifstream ifs_Y(file_Y);

    while(getline(ifs_X, name_X)){
        name_vector.push_back(name_X);
    }
    while(getline(ifs_Y, name_Y)){
        name_vector.push_back(name_Y);
    }

    int from, to, val;
    // Target type
	while(getline(ifs_WXY, str)){
		sscanf(str.c_str(), "%d %d %d", &from, &to, &val);
        struct edge_property a;
        a.label = "target";
        a.weight = val;
        property_vector.push_back(a);
        from += N;
        Wsum += val;
		edge_vector.push_back(edge(from, to));
	}
    // Attribute Type
    while(getline(ifs_WYY,str)){
		sscanf(str.c_str(), "%d %d %d", &from, &to, &val);
        struct edge_property a;
        a.label = "attribute";
        a.weight = val;
        property_vector.push_back(a);
        from += N;
        to += N;
        Wsum += val;
		edge_vector.push_back(edge(from, to));
	}
    // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
    auto tag = boost::edges_are_unsorted_multi_pass;
    // グラフのコンストラクタ
    // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
    graph g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), 5639+20);

    return g;
}

// グラフの初期化
void init_graph(graph& g)
{
    // ノードの走査
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        // (*vertex_iterator) で vertex_descriptor になる
        // g[vertex_descriptor].属性（vertex_property で定義したもの）で参照できる
        //name の入れ方name_vectorがグローバル変数
        if(*i < N ){
            g[*i].label = "target";
            g[*i].rx = 0;
            g[*i].name = name_vector[*i];
            g[*i].belongs_to_cluster = -1;
        } else {
            g[*i].label = "attribute";
            g[*i].ry = 0;
            g[*i].name = name_vector[*i];
        }
        // g[*i].previous_rank = 1.0/num_vertices(g);
        // g[*i].next_rank = 1.0/num_vertices(g);
        //g[*i].previous_rank = 0;
        //g[*i].next_rank = 0;
        g[*i].int_descriptor = static_cast<int>(*i);
    }
}

void print_detail(graph& g)
{
    // cout << "# vertices: " << num_vertices(g) << endl;
    // cout << "# edges   : " << num_edges(g) << endl;

    vertex_iterator i, j;
    // cout << " --- adjacent vertices --- " << endl;
    // for (boost::tie(i, j) = vertices(g); i!=j; i++) {
    //     cout <<  *i << " --> " << flush;
    //     // adjacent_vertices(vertex_descriptor, graph) で隣接ノードを取得できる
    //     for (auto itr = adjacent_vertices(*i, g); itr.first!=itr.second; itr.first++) {
    //         cout << *itr.first << ", " << flush;
    //     }
    //     cout << endl;
    // }


    // // out_edges(vertex_descriptor, graph) でノードの出エッジが取得できる
    // cout << " --- edges --- " << endl;
    // for (boost::tie(i, j) = vertices(g); i!=j; i++) {
    //     for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
    //         // target(edge_descriptor, graph) でエッジの先のノードを取得できる
    //         // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
    //         // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
    //         cout << *i << " --> " << target(*e.first, g) << ": " << g[*e.first].weight << " (" << g[*e.first].label << ")" << endl;
    //     }
    // }

    cout << "ranking" << endl;
    for (boost::tie(i, j) = vertices(g); *i< N ; i++) {
        //if(*i < N){
            cout << g[*i].name << "[" << g[*i].belongs_to_cluster <<"]: " << g[*i].rx << endl;
        //}else{
            //cout << *i << ":" << g[*i].ry << endl;
        //}
    }
}


