#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <chrono>
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
void ranking(graph& subgraph, int clusterNum);
void conditional_ranking(graph& g, graph& subgraph);
void get_intial_partitions(graph& g);
void write_result(vector<graph>& sub_g, string out_file);

const int iterNum = 15;
extern int xNum;
string path;
int K;
string out_file;
vector<vector<int>> cluster_label;

bool convflag = false;
int t;

int main(int argc, char* argv[])
{
    if(argc != 4){
        cout << "Error! This program needs [File Path] and [Cluster Number] [Out File]" << endl;
		cout << "Usage: " << argv[0] << "[File Path] [Cluster Number] [Out File]" << endl;
		exit(1);
	}else{
        path = argv[1];
		K = atoi(argv[2]);
        cluster_label = vector<vector<int>> (K);
        out_file = argv[3];
		if(K <= 0){
			cout << "Error: Please enter the number of clusters is 0 or more" << endl; 
			exit(0);
		}
	}

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    get_intial_partitions(g);
    cout << "< initial cluster >" << endl;
    //print_cluster(g);
    //int t;
    vector<graph> subgraph;
    for(t = 0; t < iterNum && convflag == false; t++){
        subgraph = construct_sub_graph(g);
        cout << "===== Iteration Number : " << t +1  << " =======" << endl;
        for(int clusterNum = 0; clusterNum < K; clusterNum++){
            init_graph(subgraph[clusterNum]);
            ranking(subgraph[clusterNum], clusterNum);
            conditional_ranking(g, subgraph[clusterNum]);
        }
        //print_cluster(g);
        clustering(g,subgraph);
        convflag = check_converge_cluster(g);
        //if(convflag || t == iterNum - 1)for(int clusterNum = 0; clusterNum < K; clusterNum++)print_rank_within_cluster(subgraph[clusterNum], clusterNum); 
        //print_graph_detail(g);
    }
    end = chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
    cout << " Time[micro]: " << elapsed << endl;
    cout << " Iteration Number: " << t << endl;
    write_result(subgraph, out_file);
}

bool check_converge_cluster(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); *i< xNum ; i++) {
        if(!g[*i].same_previous_cluster) return false;
    }
    cout << "converge" << endl;
    return true;
}

void conditional_ranking(graph& g, graph& subgraph){
    vertex_iterator i,j;
    double ranksum = 0;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        double tmp = 0;
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
             tmp += g[*e.first].weight * subgraph[source(*e.first, g)].ry; 
        }
        subgraph[*i].conditional_rank = tmp;
        ranksum += tmp;
    }
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        subgraph[*i].conditional_rank /= ranksum;
    }
}