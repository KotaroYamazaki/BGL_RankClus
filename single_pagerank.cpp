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
vector<double> rankSum;

void gauss_southwell(graph& g, int clusterNum);
void normalize_outedge_weight(graph& g, int clusterNum);
void normalize_outedge_weight(graph& g, int clusterNum, vector<vector<double>>& P_diff);
void normalize_xy_rank(graph& g, int clusterNum);
void normalize_global_rank(graph& g, int clusterNum);
void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);



void ranking(graph& subgraph, int clusterNum){
    if(iteration_num == 1){
        if(t == 0){
            //epsi = 1.0/(xNum + yNum);
            epsi =0;
            //cout << epsi << endl;
            normalize_outedge_weight(subgraph, clusterNum);
            single_pagerank(subgraph, clusterNum);
            pre_graph.push_back(subgraph);
        }else{
            normalize_global_rank(subgraph, clusterNum);
            normalize_outedge_weight(subgraph, clusterNum);
            init_residual(subgraph, clusterNum);
            gauss_southwell(subgraph, clusterNum);
            normalize_xy_rank(subgraph, clusterNum);
            pre_graph[clusterNum] = subgraph;
        }
    }else{
        //authority_ranking(subgraph, clusterNum);
            normalize_outedge_weight(subgraph, clusterNum);
            single_pagerank(subgraph, clusterNum);
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
            }else{
                //if(isnan(g[*i].ry)) cout << g[*i].ry << endl;
                if(RySum != 0)g[*i].ry /= RySum;
            }
        }
        rankSum[clusterNum] = (RxSum + RySum);
}

void normalize_global_rank(graph& g, int clusterNum){
    vertex_iterator i,j;
    if(g[*i].int_descriptor < xNum)g[*i].rx /= rankSum[clusterNum];
        else g[*i].ry /= rankSum[clusterNum];
//        if(isnan(g[*i].ry)) cout << g[*i].ry << endl;
}

void gauss_southwell(graph& g, int clusterNum){
    chrono::system_clock::time_point  start, end; 
    //auto msec = 0;
    queue<int> q_index;
    for(int i = 0; i < residual[clusterNum].size(); i++){
        //if(isnan(residual[clusterNum][i])) cout <<  << endl;
        if(abs(residual[clusterNum][i]) > epsi){
        //if(residual[clusterNum][i] < epsi){
            //cout << residual[clusterNum][i] << endl;
            q_index.push(i);
        }
    }
    //cout <<  "queue size : "<< q_index.size() << endl;
    for(int v = 0; v < gauss_itr; v++){
        if(q_index.empty()){
            //cout << "converged at" << v + 1 << endl;
            break;
        }
        //start = std::chrono::system_clock::now();
        unsigned long index = q_index.front();
        //cout << "index is " << index << endl;
        q_index.pop();
        vertex_descriptor max_index = vertex(index, g);
        if(isnan(max_index)){
            cout << "max_index is nan" <<endl;
            exit(0);
        }
        //cout << "res [maxindex] : "  << residual[clusterNum][max_index] << endl;
        double r_i = residual[clusterNum][max_index];
        // end = std::chrono::system_clock::now();
        // auto dur = end - start;
        // msec += chrono::duration_cast<std::chrono::microseconds>(dur).count();
        if(isnan(r_i)){
            cout << " ri: "<< r_i << endl;
            exit(1);
        }
        if(max_index < xNum){
            //cout << g[max_index].rx << endl;
            g[max_index].rx += r_i;
            // if(isnan(g[max_index].rx))exit(1);
            //cout << g[max_index].rx << endl;
        }else{
            //cout << g[max_index].ry << endl;
            if(isnan(g[max_index].ry)){
                cout << g[max_index].ry << endl;
                exit(0);
            }
            g[max_index].ry += r_i;
            //cout << g[max_index].ry << endl;
        }
        residual[clusterNum][max_index] -= r_i;
        //cout << "minus res [maxindex] : "  << residual[clusterNum][max_index] << endl;
        //double sum = 0;
        for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
           residual[clusterNum][source(*e.first, g)] +=  alpha * (g[*e.first].weight * r_i);
           //cout << target(*e.first, g) << endl;
           //cout << g[*e.first].weight << endl;
           //sum += g[*e.first].weight;
            if(abs(residual[clusterNum][target(*e.first, g)]) > epsi ) q_index.push(target(*e.first, g));
        }
        //cout << sum << endl;
        //cout << "minus2 res [maxindex] : "  << residual[clusterNum][max_index] << endl;
    }
    //cout << "Tansaku Time[micro] : " << msec << endl;
}

void init_residual(graph& g, int clusterNum){
    if(t == 1 && clusterNum == 0){
        residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
    }
        
        vertex_iterator i,j;
        for(boost::tie(i,j) = vertices(g); i != j; i++){
            //if(*i < xNum)cout << g[*i].rx << endl;
            double tmp_res = 0;
            //コスト削減のためのアドレス渡し
            graph& pre_g = pre_graph[clusterNum];
            //ランク地の引き継ぎ
            (g[*i].int_descriptor < xNum ? g[*i].rx = pre_g[*i].rx : g[*i].ry = pre_g[*i].ry);
            //if(isnan(g[*i].ry))cout << "ry is nan" << endl;

            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                vertex_descriptor v = source(*e.first, g);
                tmp_res += g[*e.first].weight * (g[v].int_descriptor < xNum ? pre_g[v].rx : pre_g[v].ry);
            }
            for (auto e = in_edges(*i, pre_g); e.first!=e.second; e.first++) {
                vertex_descriptor pre_v = source(*e.first, pre_g);
                tmp_res += pre_g[*e.first].weight * (-1)*(pre_g[pre_v].int_descriptor < xNum ? pre_g[pre_v].rx : pre_g[pre_v].ry);
            }
            //if(isnan(tmp_res))cout << "tmp_res is nan" << endl;

            residual[clusterNum][g[*i].int_descriptor] += alpha * tmp_res;
            //if(isnan(residual[clusterNum][g[*i].int_descriptor])) cout << "endl" << endl;
        }
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
                ((RxSum != 0 && g[*i].belongs_to_cluster == clusterNum) ? g[*i].rx /= RxSum : g[*i].rx = 0);
            }else{
                ((RySum != 0) ? g[*i].ry /= RySum : g[*i].ry = 0);
            }
        }
        rankSum = vector<double>(K);
        rankSum[clusterNum] = (RxSum + RySum);
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