#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;
extern vector<int> WkXY_sum;
extern int xNum;

const double alpha = 0.95;
const int rankiter = 15;

extern int K;
vector<graph> pre_graph(K);
vector<vector<double>> residual;
vector<vector<double>> pre_residual;
void gauss_southwell(graph& g, int clusterNum);

const double epsi = 0.000001;
int c = 0;


void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);

void ranking(graph& subgraph, int clusterNum){
    // if(c >= 0){
         single_pagerank(subgraph, clusterNum);
    // }else{
    //    // gauss_southwell(subgraph, clusterNum);
    // }
    //authority_ranking(subgraph, clusterNum);

    // pre_graph[clusterNum] = subgraph;
    // c++;
}

//vector<double> calc_residual(grap)
void gauss_southwell(graph& g, int clusterNum){
    vector<double> pre_rank;
    vector<double> pre_res;

    init_residual(g, clusterNum);

    for(int v = 0; v < 2000; v++){
        if(v == 0)pre_rank = init_rank(pre_graph[clusterNum]);
        //pick the largest r_i
        auto max_index = max_element(residual[clusterNum].begin(), residual[clusterNum].end());
        double r_i = residual[clusterNum][*max_index];
        if(residual[clusterNum][*max_index] < epsi)break;    
        vertex_iterator i,j;
        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(*i < xNum){
                g[*i].rx = pre_rank[*i];
                residual[clusterNum][*i] = pre_res[*i];
                if(*i == *max_index){
                    g[*i].rx += r_i;
                    residual[clusterNum][*i] -= r_i; 
                }
            }else{
                g[*i].ry = pre_rank[*i];
                if(*i == *max_index){
                    g[*i].ry += r_i;
                    residual[clusterNum][*i] -= r_i; 
                }
            }
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    if(source(*e.first, g) == *max_index)residual[clusterNum][*i] += g[*e.first].weight * r_i;
                }
        }
    }
}

vector<double> init_rank(graph& g){
    vertex_iterator i,j;
    vector<double> r;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum){
            r.push_back(g[*i].rx);
        }else{
            r.push_back(g[*i].ry);
        }
    }
    return r;
}

void init_residual(graph& g, int clusterNum){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        double tmp = 0;
        if(*i < xNum){
            g[*i].rx = pre_graph[clusterNum][*i].rx;
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    // residual 後半の項
                    tmp += (g[*e.first].weight - pre_graph[clusterNum][*e.first].weight) * pre_graph[clusterNum][source(*e.first, g)].ry;
            }
            residual[clusterNum][*i] = pre_residual[clusterNum][*i] + tmp;
        }else{
            g[*i].ry = pre_graph[clusterNum][*i].ry;
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                    // ノード *i の入エッジの重み（g[*e.first].weight）と
                    // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                    if(g[target(*e.first, g)].label == "target"){
                        tmp += (g[*e.first].weight - pre_graph[clusterNum][*e.first].weight)  * pre_graph[clusterNum][target(*e.first, g)].rx;
                    }else{
                        tmp += (g[*e.first].weight - pre_graph[clusterNum][*e.first].weight) * pre_graph[clusterNum][source(*e.first, g)].ry;
                    }
                }
            residual[clusterNum][*i] = pre_residual[clusterNum][*i] + tmp;
        }
    }
}

void single_pagerank(graph& g, int clusterNum){
    vertex_iterator i,j;

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        //出る全てのエッジに対してループを回す
        double tmp = 0;
        if(*i > xNum){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {

                    tmp += g[*e.first].weight * g[source(*e.first, g)].ry;
            }
            //residual[clusterNum][*i] = tmp;
        }else{
              for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                    // ノード *i の入エッジの重み（g[*e.first].weight）と
                    // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                    if(g[target(*e.first, g)].label == "target"){
                        tmp += (g[*e.first].weight * g[target(*e.first, g)].rx);
                    }else{
                        tmp += (g[*e.first].weight * g[source(*e.first, g)].ry);
                    }
                }
        }
    }

    //authority ranking
    double RxSum, RySum;
    for(int q = 0; q < rankiter; q++){
        RxSum = 0;
        RySum = 0;
        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            double tmp = 0;
            //出る全てのエッジに対してループを回す
            if(*i < xNum){
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
                        tmp += (g[*e.first].weight * g[target(*e.first, g)].rx);
                    }else{
                        tmp += (g[*e.first].weight * g[source(*e.first, g)].ry);
                    }
                }
                RySum += tmp;
                g[*i].ry = tmp;
            }
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(*i < xNum){
                g[*i].rx /= RxSum;
            }else{
                g[*i].ry /= RySum;
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
        if(*i < xNum){
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
            if(*i < xNum){
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
            if(*i < xNum){
                g[*i].rx /= RxSum;
            }else{
                g[*i].ry /= RySum;
            }
        }
    }
}