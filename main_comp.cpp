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
int input_seed;
string out_file;

vector<vector<int>> cluster_label;

bool convflag = false;
int iteration_num =0;

int main(int argc, char* argv[])
{
    if(argc < 3){
        cout << "Error! This program needs [File Path] and [Cluster Number] [initial seed]" << endl;
		cout << "Usage: " << argv[0] << "[File Path] [Cluster Number] [Out File]" << endl;
		exit(1);
	}else{
        path = argv[1];
		K = atoi(argv[2]);
        cluster_label = vector<vector<int>> (K);
        if(argc >= 3)input_seed = atoi(argv[3]);
        //if(argc > 3)out_file = argv[3];
		if(K <= 0){
			cout << "Error: Please enter the number of clusters is 0 or more" << endl; 
			exit(0);
		}
	}
        vector<int> time;
        for(int i = 0; i  < 2; i++){
            if(i == 0){cout << "############# <RankClus> #############" << endl;}
            else{cout << "############# Proposal #############" << endl;}
            time.push_back(do_main());
        }

	ofstream file1;
        file1.open("result_time_compare.csv",ios_base::app);
	int comp_num = 2;	
        for(int i = 0; i< comp_num;i++){
		file1 << time[i] << flush;
		if(i != comp_num - 1)file1  << "," << flush;
	}
        file1.close();

        cout << "Proposal time[micro]: " << time[1] << endl;
        cout<< "RankClus Time[micro]: " << time[0] << endl;
        cout << "Difference: " << time[0] - time[1] << endl;
        cout << "Ratio: " << 1.0*time[1]/time[0] << endl;
        cout << "NMI: " << flush;
        system("python NMI.py");
	}   

    


int do_main(){
    convflag = false;
    chrono::system_clock::time_point start, end,init_start, init_end, ranking_start, ranking_end, clustering_start, clustering_end;
    // グラフの構築
    graph g = construct_graph();

    start = chrono::system_clock::now();
    // グラフの属性値を初期化
    init_start = chrono::system_clock::now();
    init_graph(g);
    get_intial_partitions(g);
    init_end = chrono::system_clock::now();
    double init_time = std::chrono::duration_cast<std::chrono::microseconds>(init_end - init_start).count();
    cout << "initialization time[micro]: " << init_time << endl;

    vector<graph> subgraph;
    for(t = 0; t < iterNum && convflag == false; t++){
        subgraph = construct_sub_graph(g);
        cout << "===== Iteration Number : " << t + 1  << " =======" << endl;
        ranking_start = chrono::system_clock::now();
        for(int clusterNum = 0; clusterNum < K; clusterNum++){
            init_graph(subgraph[clusterNum], g);
            ranking(subgraph[clusterNum], clusterNum);
            conditional_ranking(g, subgraph[clusterNum]);
            // print_graph_detail(subgraph[clusterNum]);
        }
        ranking_end = chrono::system_clock::now();
        double ranking_time = std::chrono::duration_cast<std::chrono::microseconds>(ranking_end - ranking_start).count();
        cout << "ranking time[micro]: " << ranking_time << endl; 

        clustering_start = chrono::system_clock::now();
        clustering(g,subgraph);
        clustering_end = chrono::system_clock::now();
        double clustering_time = std::chrono::duration_cast<std::chrono::microseconds>(clustering_end - clustering_start).count();
        cout << "clustering time[micro]: " << clustering_time << endl; 
        convflag = check_converge_cluster(g);
        //print_cluster(g);
    }
    print_cluster(g);
    end = chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
    cout << " Time[micro]: " << elapsed << endl;
    cout << " Iteration Number: " << t << endl;
    cout << endl;
    //write_result(subgraph, out_file);
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
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum ; i++) {
        file << g[*i].belongs_to_cluster << flush;
        if(g[*i].int_descriptor < xNum - 1) file << "," << flush;
    }
    iteration_num++;
}

bool check_converge_cluster(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum ; i++) {
        if(!g[*i].same_previous_cluster) return false;
    }
    cout << "---------- CLUSTERS HAVE BEEN CONVERGED ---------" << endl;
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

