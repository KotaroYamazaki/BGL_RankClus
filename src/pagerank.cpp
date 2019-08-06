#include <iostream>
#include <vector>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include "gauss_southwell.hpp"
#include <chrono>

using namespace std;

extern int xNum;
extern int yNum;
extern int K;
extern int t;
extern vector<vector<int>> cluster_label;

const double alpha = 0.85;
const int rankiter = 10000;
extern double epsi;
const int gauss_start = 1;

void normalize_outedge_weight(graph& g, int clusterNum);
void pagerank_from_scratch(graph& g, int clusterNum);
void get_rank_for_rankclus(graph& g, vector<double>& r, int cluster_label);
void get_rank_for_rankclus(graph& g, int cluster_label);

// for debug
void print_cluster_with_label(graph& g);
void print_rank_within_cluster(graph& g, int clusterNum);
void print_x_p_rank(graph& g);

void ranking(graph& g, int clusterNum, int iteration_num){
    normalize_outedge_weight(g, clusterNum);
    if(iteration_num == 0){
        pagerank_from_scratch(g, clusterNum);
        get_rank_for_rankclus(g, clusterNum);
    }else{
        if(t < gauss_start){
            pagerank_from_scratch(g, clusterNum);
            get_rank_for_rankclus(g, clusterNum);
            gauss_southwell gs(epsi, alpha, clusterNum);
            gs.init(g);
        }else{
            gauss_southwell gs(epsi, alpha, clusterNum);
            gs.solve(g);
            get_rank_for_rankclus(g, clusterNum);
            gs.update_pregraph(g,clusterNum);
        }   
    }
}

void print_x_p_rank(graph& g){
    vertex_iterator i, j;
    for(boost::tie(i,j) = vertices(g); g[*i].int_descriptor < xNum; i++){
        cout << *i << ": " << g[*i].p_rank << endl;
    }
}


void normalize_outedge_weight(graph& g, int clusterNum){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                g[*e .first].weight /= out_degree(*i, g);
            }
    }
}

void pagerank_from_scratch(graph& g, int clusterNum){
    // bias value
    double b = 1.0/(xNum+yNum);
    vector<double> rank;
    rank = vector<double>(xNum + yNum, b);
    vector<double> tmp_rank = rank;
    
    vertex_iterator i,j;
    int v;
    bool conv_flag = false;
    //for(v = 0; v < rankiter; v++){
    while(!conv_flag){
	conv_flag = true;
	double change = 0;
        double ranksum = 0;

        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                tmp_rank[*i] += g[*e.first].weight * rank[source(*e.first, g)];
            }
            tmp_rank[*i] = alpha * tmp_rank[*i] + (1 - alpha) * b;
            ranksum += tmp_rank[*i];
        }

        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            tmp_rank[*i] /= ranksum;

            change = fabs(rank[*i] - tmp_rank[*i]);
            if(change > epsi)conv_flag = false;
            
            rank[*i] = tmp_rank[*i];
            g[*i].p_rank = rank[*i];
            tmp_rank[*i] = 0;
        }
        if(conv_flag)break;
    }

    //cout << "ranking converged at " << v << endl;
}

void get_rank_for_rankclus(graph& g, int cluster_label){
    //Calc rankscore for rankclus.(Convert single-graph rank score -> bi-type network graph.)
    vertex_iterator i, j;
    double rxsum = 0;
    double rysum = 0;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(cluster_label)){
            rxsum += g[*i].p_rank;
        }else{
            rysum += g[*i].p_rank;
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(cluster_label)){
            g[*i].rx = g[*i].p_rank/rxsum;
        }else{
            g[*i].ry = g[*i].p_rank/rysum;
        }
    }
}
