#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <random>
using namespace std;

graph construct_graph();
vector<graph> construct_sub_graph(graph& g);
void init_graph(graph& g);
void init_subgraph(graph& g, int clusterNum);
void print_graph_detail(graph& g);
void print_rank_within_cluster(graph& g, int clusterNum);
graph clustering(graph& g, vector<graph>& subgraph);
graph within_cluster_ranking(graph& g, int clusterNum);
void conditional_ranking(graph& g, graph& subgraph);
void get_intial_partitions(graph& g);
const int iterNum = 10;
extern int xNum;
int K = 4;

int main()
{
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    get_intial_partitions(g);
    //cout << g[boost::graph_bundle] << endl;
    vector<graph> subgraph = construct_sub_graph(g);


    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        // init_subgraph(subgraph[clusterNum], clusterNum);
        init_graph(subgraph[clusterNum]);
        within_cluster_ranking(subgraph[clusterNum], clusterNum);
        conditional_ranking(g, subgraph[clusterNum]);
        cout << "--- cluster " << clusterNum +  1 << "----" << endl;
        print_rank_within_cluster(subgraph[clusterNum], clusterNum);
    }

    clustering(g,subgraph);
    // グラフの詳細を出力
	
    // for(int t = 0; t < iterNum; t++){
    // }
}

void get_intial_partitions(graph& g){
    random_device rnd;
    vertex_iterator i, j;
    vector<bool> check_flag(K,false);
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++){
		int num = rnd()%K;
		//if(label[tmp].size() < m/K)//偏らないように調整
		//{label[tmp].push_back(i);}
		//else {i -= 1;}
        g[*i].belongs_to_cluster = num;
        check_flag[num] = true;
    }

    //check
    for(auto itr = check_flag.begin(); itr != check_flag.end(); itr++){
        if(check_flag[*itr] == false){
            get_intial_partitions(g);
        }
    }
}