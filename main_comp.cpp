#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <chrono>
#include <stdlib.h>
using namespace std;

graph construct_graph();
vector<graph> construct_sub_graph(graph& g);
void init_graph(graph& g);
void init_graph(graph& g, graph& global_g);
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

void write_result_for_NMI(graph& g);
int do_main();

const int iterNum = 15;
extern int xNum;
string path;
int K;
int t;
string out_file;

vector<vector<int>> cluster_label;

bool convflag = false;
int iteration_num =0;

int main(int argc, char* argv[])
{
    if(argc < 3){
        cout << "Error! This program needs [File Path] and [Cluster Number] [Out File]" << endl;
		cout << "Usage: " << argv[0] << "[File Path] [Cluster Number] [Out File]" << endl;
		exit(1);
	}else{
        path = argv[1];
		K = atoi(argv[2]);
        cluster_label = vector<vector<int>> (K);
        if(argc > 3)out_file = argv[3];
		if(K <= 0){
			cout << "Error: Please enter the number of clusters is 0 or more" << endl; 
			exit(0);
		}
	}
        vector<int> time;
        for(int i = 0; i  < 2; i++){
            time.push_back(do_main());
        }
        cout << "proposal time[mili] " << time[0] << endl;
        cout<< "RankClus Time[mili]: " << time[1] << endl;
        cout << "Difference " << time[1] - time[0] << endl;
        cout << "NMI: " << flush;
        system("python NMI.py");
	}   

    


int do_main(){
    convflag = false;
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // グラフの構築
    graph g = construct_graph();
    // グラフの属性値を初期化
    init_graph(g);
    get_intial_partitions(g);
    cout << "< initial cluster >" << endl;
    //print_cluster(g);

    vector<graph> subgraph;
    for(t = 0; t < iterNum && convflag == false; t++){
        subgraph = construct_sub_graph(g);
        cout << "===== Iteration Number : " << t +1  << " =======" << endl;
        for(int clusterNum = 0; clusterNum < K; clusterNum++){
            init_graph(subgraph[clusterNum], g);
            ranking(subgraph[clusterNum], clusterNum);
            conditional_ranking(g, subgraph[clusterNum]);
            // cout << "cluster: " << clusterNum << endl;
            // print_graph_detail(subgraph[clusterNum]);
        }
        //print_cluster(g);
        clustering(g,subgraph);
        convflag = check_converge_cluster(g);
       // if(convflag || t == iterNum - 1)for(int clusterNum = 0; clusterNum < K; clusterNum++)print_rank_within_cluster(subgraph[clusterNum], clusterNum); 
        //print_graph_detail(g);
    }
    end = chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
    cout << " Time[micro]: " << elapsed << endl;
    cout << " Iteration Number: " << t << endl;
    write_result(subgraph, out_file);
    write_result_for_NMI(g);
    return elapsed;
}

void write_result_for_NMI(graph& g){
    string filename;
    if(iteration_num == 0)filename = "correct.csv";
    if(iteration_num == 1)filename = "result.csv";

    fstream file;
    file.open(filename,ios::out);

    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); *i< xNum ; i++) {
        file << g[*i].belongs_to_cluster << flush;
        if(*i < xNum - 1) file << "," << flush;
    }
    iteration_num++;
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
            //if(!isnan(subgraph[source(*e.first, g)].ry))cout  << source(*e.first, g) << ": " << subgraph[source(*e.first, g)].ry << endl;
        }
        subgraph[*i].conditional_rank = tmp;
        ranksum += tmp;
    }
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        subgraph[*i].conditional_rank /= ranksum;
        //cout << subgraph[*i].conditional_rank  << endl;
    }
}
