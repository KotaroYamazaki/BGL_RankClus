#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

graph construct_graph();
vector<graph> construct_sub_graph(graph& g);
void init_graph(graph& g);
void init_subgraph(graph& g, int clusterNum);
void print_graph_detail(graph& g);
void print_rank_within_cluster(graph& g, int clusterNum);
void print_cluster(graph &g);
void clustering(graph& g, vector<graph>& subgraph);
bool check_converge_cluster(graph& g);
graph within_cluster_ranking(graph& g, int clusterNum);
void conditional_ranking(graph& g, graph& subgraph);
void get_intial_partitions(graph& g);
const int iterNum = 5;
extern int xNum;
int K = 4;
bool convflag = false;

int main()
{
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    get_intial_partitions(g);
    cout << "< initial cluster >" << endl;
    print_cluster(g);
    
    for(int t = 0; t < iterNum && convflag == false; t++){
        vector<graph> subgraph = construct_sub_graph(g);
        cout << "=====iterNum: " << t +1  << " =======" << endl;
        for(int clusterNum = 0; clusterNum < K; clusterNum++){
            init_graph(subgraph[clusterNum]);
            within_cluster_ranking(subgraph[clusterNum], clusterNum);
            conditional_ranking(g, subgraph[clusterNum]);

            //cout << "--- cluster " << clusterNum +  1 << "----" << endl;
            //print_rank_within_cluster(subgraph[clusterNum], clusterNum);
            
        }
        print_cluster(g);
        clustering(g,subgraph);
        convflag = check_converge_cluster(g);
    }

}

bool check_converge_cluster(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); *i< xNum ; i++) {
        if(!g[*i].same_previous_cluster) return false;
    }
    cout << "converge" << endl;
    return true;
}

