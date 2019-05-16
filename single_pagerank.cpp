#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <chrono>
#include <queue>

using namespace std;
extern vector<double> WkXY_sum;
extern int xNum;
extern int yNum;

const double alpha = 0.85;
const int rankiter = 15;
const int gauss_itr = 20000;
const double epsi = 0.00001;

extern int K;
vector<graph> pre_graph;
vector<vector<double>> residual;
vector<vector<vector<double>>> P;
vector<vector<vector<double>>> preP;

extern int t;
extern int iteration_num;
extern vector<vector<int>> cluster_label;

void gauss_southwell(graph& g, int clusterNum);
void normalize_outedge_weight(graph& g);
void normalize_outedge_weight(graph& g, int clusterNum, vector<vector<double>>& P_diff);
void normalize_xy_rank(graph& g, int clusterNum);
void normalize_global_rank(graph& g);
void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);

double time_init_res;
double time_gauss;

void ranking(graph& subgraph, int clusterNum){
    //clock_t start = clock();
    chrono::system_clock::time_point  start, end; 
    if(iteration_num == 1){
        if(t == 0){
            //start = std::chrono::system_clock::now();
            normalize_outedge_weight(subgraph);
            // end = std::chrono::system_clock::now();
            // auto dur = end - start;
            // auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
            // cout << "Normalize time [micro]: "<< msec << endl;
            //authority_ranking(subgraph, clusterNum);
            //start = std::chrono::system_clock::now();
            single_pagerank(subgraph, clusterNum);
            //authority_ranking(subgraph, clusterNum);
            // end = std::chrono::system_clock::now();
            // dur = end - start; 
            // msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
            // cout << "Pagerank time [micro]: "<< msec << endl;
            pre_graph.push_back(subgraph);
        }else{
            normalize_global_rank(subgraph);

            //start = std::chrono::system_clock::now();
            init_residual(subgraph, clusterNum);
            //end = std::chrono::system_clock::now();
            //dur_init_residual += end - start; 
            // auto dur = end - start; 
            // auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
            // time_init_res += msec;

            //cout <<   "init time [micro]: "<< msec << endl << endl;

            // start = std::chrono::system_clock::now();
            gauss_southwell(subgraph, clusterNum);
            // end = std::chrono::system_clock::now();
            // //dur_gauss += end - start; 
            // dur = end - start;
            // msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
            // time_gauss += msec;
            // cout <<   "gauss time [micro]: "<< msec << endl << endl;

            normalize_xy_rank(subgraph, clusterNum);
            pre_graph[clusterNum] = subgraph;
        }
        // if(clusterNum == K - 1){
        //     //auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur_init_residual).count();
        //     cout <<   "\tinit time [micro]: "<< time_init_res << endl;
        //     //auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur_gauss).count();
        //     cout <<   "\tgauss time [micro]: "<< time_gauss << endl;
        //     time_init_res = 0;
        //     time_gauss = 0;
        // }
    }else{
        authority_ranking(subgraph, clusterNum);
    }
}

void normalize_outedge_weight(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            double rowsum = 0;
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                rowsum += g[*e .first].weight;
            }
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                if(rowsum != 0)g[*e .first].weight /= rowsum;
            }
    }
}


void normalize_xy_rank(graph& g, int clusterNum){
    vertex_iterator i,j;
    double RxSum = 0;
    double RySum = 0;

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(g[*i].int_descriptor < xNum){
            if(g[*i].belongs_to_cluster == clusterNum)RxSum += g[*i].rx;
            else g[*i].rx = 0;
        }else{
            RySum += g[*i].ry;
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(g[*i].int_descriptor < xNum){
                if(RxSum != 0)g[*i].rx /= RxSum;
                //else if(g[*i].belongs_to_cluster == clusterNum)cout << "RxSum == 0: "<< g[*i].rx  << endl;
                //cout<< g[*i].rx << endl;
            }else{
                if(RySum != 0)g[*i].ry /= RySum;
            }
        }
}

void normalize_global_rank(graph& g){
    vertex_iterator i,j;
    double ranksum = 0;
    double rank;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {   
        if(g[*i].int_descriptor < xNum)rank = g[*i].rx;
        else rank = g[*i].ry;
        ranksum += rank;
    }
    if(g[*i].int_descriptor < xNum)g[*i].rx /= ranksum;
        else g[*i].ry /= ranksum;

}

void gauss_southwell(graph& g, int clusterNum){
    chrono::system_clock::time_point  start, end; 
    //auto msec = 0;
    queue<int> q_index;
    for(int i = 0; i < residual[clusterNum].size(); i++){
        if(residual[clusterNum][i] < epsi){
            q_index.push(i);
        }
    }
    for(int v = 0; v < gauss_itr; v++){
        if(q_index.empty()){
            //cout << "converged at" << v + 1 << endl;
            break;
        }
        //start = std::chrono::system_clock::now();
        unsigned long index = q_index.front();
        q_index.pop();
        vertex_descriptor max_index = vertex(index, g);
        double r_i = residual[clusterNum][max_index];
        // end = std::chrono::system_clock::now();
        // auto dur = end - start;
        // msec += chrono::duration_cast<std::chrono::microseconds>(dur).count();

        if(max_index < xNum){
            //cout << g[max_index].rx << endl;
            g[max_index].rx += r_i;
            // if(isnan(g[max_index].rx))exit(1);
            //cout << g[max_index].rx << endl;
        }else{
            g[max_index].ry += r_i;
        }
        residual[clusterNum][max_index] -= r_i;
        //double sum = 0;
        for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
           residual[clusterNum][target(*e.first, g)] +=  alpha * (g[*e.first].weight * r_i);
           //sum += g[*e.first].weight;
            if(residual[clusterNum][target(*e.first, g)] > epsi ) q_index.push(target(*e.first, g));
        }
        //cout << sum << endl;
    }
    //cout << "Tansaku Time[micro] : " << msec << endl;
}

void init_residual(graph& g, int clusterNum){
    chrono::system_clock::time_point  start, end; 
    
    if(t == 1 && clusterNum == 0){
        residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
        //P_diff = vector<vector<double>>(num_vertices(g), vector<double>(num_vertices(g), 0));
    }
        //start = std::chrono::system_clock::now();
        // //normalize_outedge_weight(g, clusterNum, P_diff);
        normalize_outedge_weight(g);
        //end = std::chrono::system_clock::now();
        // auto dur = end - start; 
        // auto msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
        // cout << "   normalize time: " << msec << endl;

        vertex_iterator i,j,m,n;
        for(boost::tie(i,j) = vertices(g); i != j; i++){
            //if(*i < xNum)cout << g[*i].rx << endl;
            double tmp_res = 0;
            //コスト削減のためのアドレス渡し
            graph& pre_g = pre_graph[clusterNum];
            //ランク地の引き継ぎ
            (g[*i].int_descriptor < xNum ? g[*i].rx = pre_g[*i].rx : g[*i].ry = pre_g[*i].ry);
            //if(g[*i].belongs_to_cluster == clusterNum)cout << pre_g[*i].rx << endl;

            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                vertex_descriptor v = source(*e.first, g);
                tmp_res += g[*e.first].weight * (g[v].int_descriptor < xNum ? pre_g[v].rx : pre_g[v].ry);
            }
            for (auto e = in_edges(*i, pre_g); e.first!=e.second; e.first++) {
                vertex_descriptor pre_v = source(*e.first, pre_g);
                tmp_res += pre_g[*e.first].weight * (pre_g[pre_v].int_descriptor < xNum ? pre_g[pre_v].rx : pre_g[pre_v].ry);
            }

            residual[clusterNum][g[*i].int_descriptor] += alpha * tmp_res;
        }

        // end = std::chrono::system_clock::now();
        // dur = end - start; 
        // msec = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
        // cout << "   Keisan time: " << msec << endl;
}

void single_pagerank(graph& g, int clusterNum){
    vertex_iterator i,j;
    double RxSum = 0;
    double RySum = 0;
    for(int v = 0; v < rankiter; v++){
        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
                if(v == 0){
                    if(g[*i].int_descriptor < xNum ){
                        g[*i].rx = 1.0/cluster_label[clusterNum].size();
                        if(g[*i].belongs_to_cluster != clusterNum)g[*i].rx = -1;
                    }else{
                        g[*i].ry = 1.0/yNum;
                    }
                }

            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                if(g[*i].int_descriptor < xNum){
                    g[*i].rx =  alpha * g[*e.first].weight * g[*i].rx+ (1.0 - alpha)* (1.0/(cluster_label[clusterNum].size()));
                    RxSum += g[*i].rx;
                    
                }else{
                    g[*i].ry = alpha * g[*e.first].weight;
                    RySum += g[*i].ry;
                }
            }
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(g[*i].int_descriptor < xNum){
                //if(RxSum != 0)g[*i].rx /= RxSum;
                ((RxSum != 0 && g[*i].belongs_to_cluster == clusterNum) ? g[*i].rx /= RxSum : g[*i].rx = 0);
                //if(g[*i].belongs_to_cluster == clusterNum)cout << *i << ": "<< g[*i].rx << ":" << g[*i].belongs_to_cluster <<endl;
            }else{
                //if(RySum != 0) g[*i].ry /= RySum;
                ((RySum != 0) ? g[*i].ry /= RySum : g[*i].ry = 0);
                //cout << g[*i].ry << endl;
            }
        }
}

void authority_ranking(graph& g, int clusterNum){
    vertex_iterator i, j;
    //iterator を用いて最初のiterから最後まで回す

    //simple ranking
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        double tmp = 0;
        //出る全てのエッジに対してループを回す
        if(g[*i].int_descriptor < xNum){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                // ノード *i の入エッジの重み（g[*e.first].weight）と
                // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                tmp += g[*e.first].weight;
            }
            g[*i].rx = tmp/WkXY_sum[clusterNum];
        }else{
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                // ノード *i の入エッジの重み（g[*e.first].weight）と
                // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                if(g[target(*e.first, g)].label == "target")tmp += g[*e.first].weight;
            }
            g[*i].ry = tmp/WkXY_sum[clusterNum];
        }
    }

    //authority ranking
    for(int q = 0; q < rankiter; q++){
        double RxSum = 0;
        double RySum = 0;
        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            double tmp = 0;
            //出る全てのエッジに対してループを回す
            if(g[*i].int_descriptor < xNum){
                for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    // ノード *i の入エッジの重み（g[*e.first].weight）と
                    // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                    tmp += g[*e.first].weight * g[source(*e.first, g)].ry;
                }
                RxSum += tmp;
                g[*i].rx = tmp;
            }else{
                for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                    // ノード *i の入エッジの重み（g[*e.first].weight）と
                    // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                    if(g[target(*e.first, g)].label == "target"){
                        tmp += alpha*(g[*e.first].weight * g[target(*e.first, g)].rx);
                    }else{
                        tmp += (1-alpha)*(g[*e.first].weight * g[source(*e.first, g)].ry);
                    }
                }
                RySum += tmp;
                g[*i].ry = tmp;
            }
        }

        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(g[*i].int_descriptor < xNum){
                g[*i].rx /= RxSum;
            }else{
                g[*i].ry /= RySum;
            }
        }
    }
}