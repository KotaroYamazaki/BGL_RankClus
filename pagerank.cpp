#include <iostream>
#include <vector>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <chrono>
#include <queue>

using namespace std;

extern int xNum;
extern int yNum;
extern int K;
extern int t;
extern int iteration_num;
extern vector<vector<int>> cluster_label;
extern vector<vector<double>> row_sum_vec;

const double alpha = 0.85;
const int rankiter = 30;
const int gauss_itr = 20000;
double epsi;

vector<graph> pre_graph;
vector<vector<double>> residual;

void normalize_outedge_weight(graph& g, int clusterNum);
void pagerank(graph& g, int clusterNum);
void get_rank_for_rankclus(graph& g, vector<double> r, int cluster_label, double rxsum, double rysum);


void ranking(graph& g, int clusterNum){
    if(iteration_num == 0){
        normalize_outedge_weight(g, clusterNum);
        pagerank(g, clusterNum);
    }
}

void normalize_outedge_weight(graph& g, int clusterNum){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                if(row_sum_vec[clusterNum][*i] != 0)g[*e .first].weight /= row_sum_vec[clusterNum][*i];
            }
    }
}

void pagerank(graph& g, int clusterNum){
    double b = 1.0/(xNum+yNum);
    vector<double> rank;
    rank = vector<double>(xNum + yNum, b);
    vector<double> tmp_rank = rank;
    
    vertex_iterator i,j;

    // 正規化はいらない
    double RxSum = 0;
    double RySum = 0;
    for(int v = 0; v < rankiter; v++){
        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                tmp_rank[*i] += g[*e.first].weight * rank[source(*e.first, g)];
            }
            tmp_rank[*i] = alpha * tmp_rank[*i] + (1 - alpha) * b;
            rank[*i] = tmp_rank[*i];
            tmp_rank[*i] = 0;

            if(v == rankiter - 1){
                (*i < xNum) ? RxSum += rank[*i] : RySum += rank[*i];
            }
        }
    }
    get_rank_for_rankclus(g, rank, clusterNum, RxSum, RySum);
}

void get_rank_for_rankclus(graph& g, vector<double> r, int cluster_label, double rxsum, double rysum){
    //Calc rankscore for rankclus.(Convert single-graph rank score -> bi-type network graph.)
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(cluster_label)){
            g[*i].rx = r[*i]/rxsum;
        }else{
            g[*i].ry = r[*i]/rysum;
        }
        g[*i].p_rank = r[*i];
    }
}