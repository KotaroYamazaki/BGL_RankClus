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
const int gauss_start = 0;

vector<graph> pre_graph;
vector<vector<double>> residual;

void normalize_outedge_weight(graph& g, int clusterNum);
void pagerank_from_scratch(graph& g, int clusterNum);
void get_rank_for_rankclus(graph& g, vector<double>& r, int cluster_label);
void get_rank_for_rankclus(graph& g, int cluster_label);
void calc_initial_residual(graph& g);
pair<queue<int>, vector<bool>> calc_tracking_residual(graph& g, int clusterNum);
void gauss_southwell(graph& g, int clusterNum, queue<int>& q, vector<bool>& occupied_flag);

// for debug
void print_cluster_with_label(graph& g);
void print_rank_within_cluster(graph& g, int clusterNum);
void print_x_p_rank(graph& g);

class PreGraph_Structure{
    vector<vector<int>> in_node_set;
    vector<double> p_rank;
};

void ranking_with_time(graph &g, int clusterNum){
    auto start_r = std::chrono::system_clock::now();
    pair<queue<int>, vector<bool>> p = calc_tracking_residual(g, clusterNum);
    auto end_r = std::chrono::system_clock::now();
    auto dur_r = end_r - start_r;        // 要した時間を計算
    auto msec_r = std::chrono::duration_cast<std::chrono::microseconds>(dur_r).count();
    std::cout << "cluster = "<<clusterNum << " residu time [micro] : "<< msec_r << "\n";

    auto start_g = std::chrono::system_clock::now(); 
    gauss_southwell(g, clusterNum, p.first, p.second);
    auto end_g = std::chrono::system_clock::now();
    auto dur = end_g - start_g;        // 要した時間を計算
    auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
    std::cout << "cluster = "<<clusterNum << " gauss time[micro] : "<< msec << endl;

    auto start_get = std::chrono::system_clock::now(); 
    get_rank_for_rankclus(g, clusterNum);
    auto end_get = std::chrono::system_clock::now();
    auto dur_get = end_get - start_get;        // 要した時間を計算
    auto msec_get = std::chrono::duration_cast<std::chrono::microseconds>(dur_get).count();
    std::cout << "cluster = "<<clusterNum << " get rank time[micro] : "<< msec_get << endl;

    auto start_p = std::chrono::system_clock::now();
    pre_graph[clusterNum] = g;
    auto end_p = std::chrono::system_clock::now();
    auto dur_p = end_p - start_p;        // 要した時間を計算
    auto msec_p = std::chrono::duration_cast<std::chrono::microseconds>(dur_p).count();
    std::cout << "cluster = "<<clusterNum << " pregraohp time[micro] : "<< msec_p << endl;
}


void ranking(graph& g, int clusterNum){
    normalize_outedge_weight(g, clusterNum);
    //epsi = 0.0001/(xNum+yNum);
    epsi = 1e-9;
    //epsi = 1e-11;
    if(iteration_num == 0){
        pagerank_from_scratch(g, clusterNum);
        get_rank_for_rankclus(g, clusterNum);
    }else{
        if(t == gauss_start){
            pagerank_from_scratch(g, clusterNum);
            get_rank_for_rankclus(g, clusterNum);

            //auto start_c = std::chrono::system_clock::now();
            if(t == gauss_start)calc_initial_residual(g);
            // auto end_c = std::chrono::system_clock::now();
            // auto dur_c = end_c - start_c;        // 要した時間を計算
            // auto msec_c = std::chrono::duration_cast<std::chrono::microseconds>(dur_c).count();
            // std::cout << "cluster = "<<clusterNum << "  calc res time [micro] : "<< msec_c << "\n";


            //auto start_r = std::chrono::system_clock::now();
            if(t == 0)pre_graph.push_back(g);
            else pre_graph[clusterNum] = g;
            // auto end_r = std::chrono::system_clock::now();
            // auto dur_r = end_r - start_r;        // 要した時間を計算
            // auto msec_r = std::chrono::duration_cast<std::chrono::microseconds>(dur_r).count();
            // std::cout << "cluster = "<<clusterNum << "  graph push time [micro] : "<< msec_r << "\n";
            
        }else{
            pair<queue<int>, vector<bool>> p = calc_tracking_residual(g, clusterNum);
            gauss_southwell(g, clusterNum, p.first, p.second);
            get_rank_for_rankclus(g, clusterNum);
            pre_graph[clusterNum] = g;
            //ranking_with_time(g, clusterNum);
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
        // for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
        //     unsigned long index = source(*e.first, g);
        //     res[index] +=  alpha * (g[*e.first].weight * r_i);
        //     if(fabs(res[index]) > epsi && !occupied_flag[index]){
        //         q.push(index);
        //         occupied_flag[index] = true;
        //     }
        // }
        for (auto e = out_edges(max_index, g); e.first!=e.second; e.first++) {
            unsigned long index = target(*e.first, g);
            res[index] +=  alpha * (g[*e.first].weight * r_i);
            if(fabs(res[index]) > epsi && !occupied_flag[index]){
                q.push(index);
                occupied_flag[index] = true;
            }
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
        //cout << tmp << endl;
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
        double Pt_times_x = 0;
        double Pt_mius_1_times_x = 0;
        //calc P(t)*x(t-1)
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            v = source(*e.first, g);
            Pt_times_x += g[*e.first].weight * pre_g[v].p_rank;
        }
        //calc P(t-1)*x(t-1)
        for (auto e = in_edges(*i, pre_g); e.first!=e.second; e.first++) {
            v = source(*e.first, pre_g);
            Pt_mius_1_times_x += pre_g[*e.first].weight * pre_g[v].p_rank;
        }
        
        residual[clusterNum][*i] += alpha * (Pt_times_x - Pt_mius_1_times_x);
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
    for(v = 0; v < rankiter; v++){
        double change = 0;
        bool conv_flag = true;
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