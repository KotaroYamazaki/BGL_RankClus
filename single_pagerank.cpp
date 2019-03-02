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

extern vector<int> WkXY_sum;
extern int K;
vector<graph> pre_graph;
vector<vector<double>> residual;
vector<vector<double>> pre_residual;
void gauss_southwell(graph& g, int clusterNum);

const double epsi = 0.000001;

extern int t;


void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);

void ranking(graph& subgraph, int clusterNum){
    if(t == 0){
        single_pagerank(subgraph, clusterNum);
        pre_graph.push_back(subgraph);
    }else{
        // subgraph = normalize_weight(subgraph, clusterNum);
        gauss_southwell(subgraph, clusterNum);
        pre_graph[clusterNum] = subgraph;
    }
}

//vector<double> calc_residual(grap)
void gauss_southwell(graph& g, int clusterNum){
    init_residual(g, clusterNum);

    for(int v = 0; v < 20000; v++){
        auto max_itr = max_element(residual[clusterNum].begin(), residual[clusterNum].end());
        //vertex_iterator m;
        //m = max_element(residual[clusterNum].begin(), residual[clusterNum].end());
        //*m = distance(residual[clusterNum].begin(), max_itr);
        unsigned long max_index = distance(residual[clusterNum].begin(), max_itr);
        double r_i = residual[clusterNum][max_index];

        if(r_i < epsi){
            cout << "converge at " << v << endl;
            break;
        }

        if(max_index < xNum){
            g[max_index].rx += r_i;
        }else{
            g[max_index].ry += r_i;
        }
        residual[clusterNum][max_index] -= r_i; 

        vertex_iterator i,j;
                for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
                    residual[clusterNum][max_index] += g[*e.first].weight * r_i;
                }
        //     }
        // }

        // pre_res = residual[clusterNum];

        // for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        //     if(*i < xNum){
        //         g[*i].rx /= RxSum;
        //     }else{
        //         g[*i].ry /= RySum;
        //     }
        // }

    }
}


bool is_reachable(vertex_descriptor s, vertex_descriptor t, graph& g ){
    for (auto e = in_edges(s, g); e.first!=e.second; e.first++) {
        if(source(*e.first, g) == t)return true;
    }
    return false;
}


void init_residual(graph& g, int clusterNum){
    vertex_iterator i,j;
    if(t == 1){
        residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
        pre_residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
    }
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        double tmp = 0;
        if(*i < xNum){
            g[*i].rx = pre_graph[clusterNum][*i].rx;
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    // residual 後半の項
                    //double pre_weight;
                    //std::vector<boost::default_color_type> color(vertices(g), boost::white_color);
                    if(!is_reachable(*i, source(*e.first, g), pre_graph[clusterNum])){
                        tmp += (g[*e.first].weight - 0)  * pre_graph[clusterNum][source(*e.first, g)].ry;
                    }
            }
            residual[clusterNum][*i] = pre_residual[clusterNum][*i] + tmp;
        }else{
            g[*i].ry = pre_graph[clusterNum][*i].ry;
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    // ノード *i の入エッジの重み（g[*e.first].weight）と
                    // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
                    if(g[source(*e.first, g)].label == "target"){
                        if(!is_reachable(*i, source(*e.first, g), pre_graph[clusterNum]))tmp += (g[*e.first].weight - 0)  * pre_graph[clusterNum][target(*e.first, g)].rx; //tmp += (g[*e.first].weight - pre_graph[clusterNum][*e.first].weight)  * pre_graph[clusterNum][target(*e.first, g)].rx;
                    }else{
                        if(!is_reachable(*i, source(*e.first, g), pre_graph[clusterNum]))tmp += (g[*e.first].weight - 0) * pre_graph[clusterNum][target(*e.first, g)].ry;
                    }
                }
            residual[clusterNum][*i] = pre_residual[clusterNum][*i] + tmp;
        }
    }
}

void single_pagerank(graph& g, int clusterNum){
    vertex_iterator i,j;

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