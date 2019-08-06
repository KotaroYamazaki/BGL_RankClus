#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <stdlib.h>
#include "graph.hpp"
#include "clustering.hpp"
#include "pagerank.hpp"
using namespace std;

void ranking(graph& subgraph, int clusterNum, int iteration_num);

const int iterNum = 15;
extern int xNum;
extern double WXY_sum;
extern vector<double> WkXY_sum;
double epsi;
string path;
int K;
int t;
int input_seed;
string out_file;

vector<vector<int>> cluster_label;

bool convflag = false;

int main(int argc, char* argv[])
{
    if(argc < 4){
        cout << "Error! This program needs [-r/-d] [File Path] and [Cluster Number] [initial seed] [epsilon^<>]" << endl;
        cout << "Usage: " << argv[0] << "[File Path] [Cluster Number] [initial seed] [epsilon^<>]" << endl;
        exit(1);
	}else{
        path = argv[1];
		K = atoi(argv[2]);
        if(argc >= 3)input_seed = atoi(argv[3]);
        if(argc > 3)epsi = atoi(argv[4]);
        epsi = pow(10,-epsi);
        cout <<"Cluster Number: " << K << endl;
        cout << "Input seed: " << input_seed << endl;
        cout << "Epsilon: " << epsi << endl;

		if(K <= 0){
			cout << "Error: Please enter the number of clusters is 0 or more" << endl; 
			exit(0);
		}
	}

        vector<int> time;
        for(int i = 0; i  < 2; i++){
            double r_time = 0;
            double c_time = 0;

            if(i == 0){cout << "############# <RankClus> #############" << endl;}
            else{cout << "############# Proposal #############" << endl;}

            cluster_label = vector<vector<int>> (K);
            convflag = false;
            chrono::system_clock::time_point start, end,init_start, init_end, ranking_start, ranking_end, clustering_start, clustering_end;
            // グラフの構築
            graph g = construct_graph();

            start = chrono::system_clock::now();
            // グラフの属性値を初期化
            init_start = chrono::system_clock::now();
            cout << "init graph." << endl;
            init_graph(g);
            cout << "get initial partition." << endl;
            get_intial_partitions(g);
            init_end = chrono::system_clock::now();
            double init_time = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start).count();
            cout << "initialization time[milli]: " << init_time << endl;

            vector<graph> subgraph;
            for(t = 0; t < iterNum && convflag == false; t++){
                subgraph = construct_sub_graph(g);
                cout << "===== Iteration Number : " << t + 1  << " =======" << endl;
                ranking_start = chrono::system_clock::now();
                for(int clusterNum = 0; clusterNum < K; clusterNum++){
                    init_graph(subgraph[clusterNum], g);
                    pagerank pr(xNum, t, clusterNum, epsi);
                    pr.solve(subgraph[clusterNum], i);
                    conditional_ranking(g, subgraph[clusterNum]);
                }
                ranking_end = chrono::system_clock::now();
                double ranking_time = std::chrono::duration_cast<std::chrono::milliseconds>(ranking_end - ranking_start).count();
                r_time += ranking_time;
                cout << "Ranking time[milli]: " << ranking_time << endl;

                clustering_start = chrono::system_clock::now();
                clustering cl(WkXY_sum, K, WXY_sum, xNum, cluster_label);
                cluster_label = cl.update_cluster_label(g, subgraph);
                clustering_end = chrono::system_clock::now();
                double clustering_time = std::chrono::duration_cast<std::chrono::milliseconds>(clustering_end - clustering_start).count();
                c_time += clustering_time;
                cout << "Clustering time[milli]: " << clustering_time << endl;
                convflag = check_converge_cluster(g);
                write_result_in_each_iteration(g, i, t);
            }
            end = chrono::system_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
            cout << " Time[milli]: " << elapsed << endl;
            //double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
            //double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
            cout << " Iteration Number: " << t << endl;
            cout << endl;
            time.push_back(r_time);
            time.push_back(c_time);
            time.push_back(elapsed);
            write_result_for_NMI(g,i);
        }
        write_result_to_csv(time);

        cout<< "RankClus Time[milli]: " << time[0] << endl;
        cout << "Proposal time[milli]: " << time[1] << endl;
        cout << "Difference: " << time[0] - time[1] << endl;
        //cout << "Ratio: " << 1.0*time[1]/time[0] << endl;
        cout << "NMI: " << flush;
        system("python src/NMI.py");
        system("python src/NMI_in_each_itr.py");
	}



