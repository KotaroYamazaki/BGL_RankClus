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
const int rankiter = 10000;
double epsi;

vector<graph> pre_graph;
vector<vector<double>> residual;

void normalize_outedge_weight(graph& g, int clusterNum);
void pagerank_from_scratch(graph& g, int clusterNum);
void get_rank_for_rankclus(graph& g, vector<double> r, int cluster_label);
void get_rank_for_rankclus(graph& g, int cluster_label);
void calc_initial_residual(graph& g);
pair<queue<int>, vector<bool>> calc_tracking_residual(graph& g, int clusterNum);
void gauss_southwell(graph& g, int clusterNum, queue<int>& q, vector<bool>& occupied_flag);

// for debug
void print_cluster_with_label(graph& g);
void print_rank_within_cluster(graph& g, int clusterNum);
void print_x_p_rank(graph& g);

void ranking(graph& g, int clusterNum){
    normalize_outedge_weight(g, clusterNum);
    //epsi = 0.1/(xNum+yNum);
    epsi = 0;
    if(iteration_num == 0){
        pagerank_from_scratch(g, clusterNum);
        get_rank_for_rankclus(g, clusterNum);
    }else{
        if(t == 0){
            pagerank_from_scratch(g, clusterNum);
            get_rank_for_rankclus(g, clusterNum);
            calc_initial_residual(g);
            pre_graph.push_back(g);
        }else{
            pair<queue<int>, vector<bool>> p = calc_tracking_residual(g, clusterNum);
            gauss_southwell(g, clusterNum, p.first, p.second);
            get_rank_for_rankclus(g, clusterNum);
            pre_graph[clusterNum] = g;
        }
        
    }
    //print_rank_within_cluster(g, clusterNum);
    //cout << "< cluster = " << clusterNum << " > "<< endl;
    //print_x_p_rank( g);
}

void gauss_southwell(graph& g, int clusterNum, queue<int>& q, vector<bool>& occupied_flag){
    vector<double>& res = residual[clusterNum];

    while(!q.empty()){
        unsigned long index = q.front();
        vertex_descriptor max_index = vertex(index, g);
        q.pop();
        occupied_flag[max_index] = false;
        
        double r_i = res[max_index];
        //x^(v) = x^(v-1) + r_i^(v-1)*e_i
        g[max_index].p_rank += r_i;

        //r^(v) = r^(v-1) - r_i^(v-1)*e_i
        res[max_index] -= r_i;

        //r^(v) = r^(v-1) + a*r_i^(v-1)*P*e_i
        for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
            unsigned long index = source(*e.first, g);
            res[index] +=  alpha * (g[*e.first].weight * r_i);
            if(occupied_flag[index] && fabs(res[index]) > epsi) q.push(index);
        }
    }
}

void calc_initial_residual(graph& g){
    vector<double> tmp_res;
    double b = 1.0/(xNum+ yNum);
    vertex_iterator i,j;
    for(boost::tie(i,j) = vertices(g); i != j; i++){
        double tmp = (1-alpha)*b;
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            tmp += alpha * g[*e.first].weight * g[source(*e.first, g)].p_rank;
        }
        tmp -= g[*i].p_rank;
        tmp_res.push_back(tmp);
    }
    residual.push_back(tmp_res);
}

void print_x_p_rank(graph& g){
    vertex_iterator i, j;
    for(boost::tie(i,j) = vertices(g); g[*i].int_descriptor < xNum; i++){
        cout << *i << ": " << g[*i].p_rank << endl;
    }
}

pair<queue<int>, vector<bool>> calc_tracking_residual(graph& g,int clusterNum){
    queue<int> q;
    vector<bool> occupied_flag = vector<bool>(xNum+yNum);
    vertex_iterator i, j;
    graph& pre_g = pre_graph[clusterNum];
    for(boost::tie(i,j) = vertices(g); i != j; i++){
        g[*i].p_rank = pre_g[*i].p_rank;
        vertex_descriptor v;
        double tmp_res = 0;
        //calc P(t)*x(t-1)
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            v = source(*e.first, g);
            tmp_res += g[*e.first].weight * g[*i].p_rank;
        }
        //calc P(t-1)*x(t-1)
        for (auto e = in_edges(*i, pre_g); e.first!=e.second; e.first++) {
            v = source(*e.first, pre_g);
            tmp_res -= pre_g[*e.first].weight * pre_g[*i].p_rank;
        }
        residual[clusterNum][*i] += alpha * tmp_res;
        if(fabs(residual[clusterNum][*i]) > epsi){
            q.push(*i);
            occupied_flag[*i] = true;
        }
    }
    return make_pair(q, occupied_flag);
}



void normalize_outedge_weight(graph& g, int clusterNum){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                if(row_sum_vec[clusterNum][*i] != 0)g[*e .first].weight /= row_sum_vec[clusterNum][*i];
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
    for(v = 0; v < rankiter; v++){
        double change = 0;
        bool conv_flag = true;
        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                tmp_rank[*i] += g[*e.first].weight * rank[source(*e.first, g)];
            }
            tmp_rank[*i] = alpha * tmp_rank[*i] + (1 - alpha) * b;
            change = fabs(rank[*i] - tmp_rank[*i]);
            if(change > epsi)conv_flag = false;
            rank[*i] = tmp_rank[*i];
            tmp_rank[*i] = 0;
        }
        if(conv_flag)break;
    }
    for(boost::tie(i,j) = vertices(g); i !=j; i++){
        //if(*i < xNum)cout << *i << ": " << rank[*i] << "<- " << tmp_rank[*i] << endl;  
        g[*i].p_rank = rank[*i];
    }
    cout << "ranking converged at " << v << endl;
}

void get_rank_for_rankclus(graph& g, vector<double> r, int cluster_label){
    //Calc rankscore for rankclus.(Convert single-graph rank score -> bi-type network graph.)
    vertex_iterator i, j;
    double rxsum = 0;
    double rysum = 0;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(cluster_label)){
            rxsum += r[*i];
        }else{
            rysum += r[*i];
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(cluster_label)){
            g[*i].rx = r[*i]/rxsum;
        }else{
            g[*i].ry = r[*i]/rysum;
        }
        g[*i].p_rank = r[*i];
    }
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