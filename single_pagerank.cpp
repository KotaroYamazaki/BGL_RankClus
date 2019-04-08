#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <time.h>
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
graph normalize_weight(graph& g);
void normalize_rank(graph& g, int clusterNum);

const double epsi = 0.000001;

extern int t;
extern int iteration_num;


void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);

void ranking(graph& subgraph, int clusterNum){
    //clock_t start = clock();
    if(iteration_num == 1){
        if(t == 0){
            authority_ranking(subgraph, clusterNum);
            pre_graph.push_back(subgraph);
        }else{
            init_residual(subgraph, clusterNum);
            gauss_southwell(subgraph, clusterNum);
            normalize_rank(subgraph, clusterNum);
            pre_graph[clusterNum] = subgraph;
        }
    }else{
        authority_ranking(subgraph, clusterNum);
    }
}



graph normalize_weight(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        //O(n*E*E)
        if(g[*i].int_descriptor < xNum){
            //g[*i].rx = pre_graph[clusterNum][*i].rx;
            int rowsum = 0;
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                rowsum += g[*e .first].weight;
            }
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                g[*e .first].weight /= rowsum;
            }
        }
    }
    return g;
}

void normalize_rank(graph& g, int clusterNum){
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
                g[*i].rx /= RxSum;
                // cout<< g[*i].rx << endl;
            }else{
                g[*i].ry /= RySum;
            }
        }

}

void gauss_southwell(graph& g, int clusterNum){
    normalize_weight(g);
    for(int v = 0; v < 20000; v++){
        auto max_itr = max_element(residual[clusterNum].begin(), residual[clusterNum].end());
        unsigned long max_index = distance(residual[clusterNum].begin(), max_itr);
        double r_i = residual[clusterNum][max_index];
        //cout << "r_i:  " <<  r_i << endl;

        if(r_i < epsi){
            // cout << "r_i: " << r_i << endl;
            // cout << "converge at " << v << endl;
            break;
        }

        if(max_index < xNum){
            g[max_index].rx += r_i;
            //cout<< g[max_index].rx << endl;
        }else{
            g[max_index].ry += r_i;
        }
        residual[clusterNum][max_index] -= r_i; 

        for (auto e = in_edges(max_index, g); e.first!=e.second; e.first++) {
            residual[clusterNum][max_index] += g[*e.first].weight * r_i;
        }
    }
}

void init_residual(graph& g, int clusterNum){
    vertex_iterator i,j;
    if(t >= 1){
        residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
        pre_residual = residual;
    }
    residual = pre_residual;
    //O(n*E)
    for (boost::tie(i, j) = vertices(g); i != j; i++) {
        if(g[*i].int_descriptor < xNum){
            if(isnan(pre_graph[clusterNum][*i].rx)) cout << "pregraph rx is nan" <<endl;
            g[*i].rx = pre_graph[clusterNum][*i].rx;
            double tmp = 0;
            // クラスタが変化している場合（変化していない場合はそのまま前回計算した残差を用いる）
            if(pre_graph[clusterNum][*i].belongs_to_cluster != g[*i].belongs_to_cluster){
                // 前は所属していたというノード
                if(pre_graph[clusterNum][*i].belongs_to_cluster == clusterNum){
                    for (auto e = in_edges(*i, pre_graph[clusterNum]); e.first!=e.second; e.first++) {
                        tmp += (0 - pre_graph[clusterNum][*e.first].weight) * pre_graph[clusterNum][source(*e.first, pre_graph[clusterNum])].ry;
                        residual[clusterNum][source(*e.first, pre_graph[clusterNum])] += (0 - pre_graph[clusterNum][*e.first].weight)*pre_graph[clusterNum][*i].rx;
                    }
                    residual[clusterNum][g[*i].int_descriptor] += tmp;
                }else{
                    //前は所属していなかったが今は所属しているノード
                    for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                        tmp += (g[*e.first].weight - 0) * pre_graph[clusterNum][source(*e.first, g)].ry;
                        residual[clusterNum][source(*e.first, g)] += (0 - pre_graph[clusterNum][*e.first].weight)*pre_graph[clusterNum][*i].rx;
                    }
                    residual[clusterNum][g[*i].int_descriptor] += tmp;
                }
            }
        }else{
            g[*i].ry = pre_graph[clusterNum][*i].ry;
        }
        //residual[clusterNum][g[*i].int_descriptor] = pre_residual[clusterNum][g[*i].int_descriptor] + tmp;
        

        // if(g[*i].int_descriptor < xNum){
        //     g[*i].rx = pre_graph[clusterNum][*i].rx;
        //     for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
        //             // residual 後半の項
        //             // double pre_weight;
        //             if(!boost::edge(*i, source(*e.first, g), pre_graph[clusterNum]).second){
        //                 tmp += (g[*e.first].weight - 0)  * pre_graph[clusterNum][source(*e.first, g)].ry;
        //             }
        //     }
        //     //if(tmp > 0 )cout << tmp << endl;
        //     residual[clusterNum][g[*i].int_descriptor] = pre_residual[clusterNum][g[*i].int_descriptor] + tmp;
        // }else{
        //     g[*i].ry = pre_graph[clusterNum][*i].ry;
        //     for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
        //             // ノード *i の入エッジの重み（g[*e.first].weight）と
        //             // そのエッジの元ノード（source(*e.first, g) のランク値（g[source(*e.first, g)].previous_rank）をかける
        //             if(g[source(*e.first, g)].label == "target" && !boost::edge(*i, source(*e.first, g), pre_graph[clusterNum]).second){
        //                 tmp += (g[*e.first].weight - 0)  * pre_graph[clusterNum][target(*e.first, g)].rx;
        //             }else{
        //                 tmp += (g[*e.first].weight - 0) * pre_graph[clusterNum][target(*e.first, g)].ry;
        //             }
        //         }
        //     residual[clusterNum][g[*i].int_descriptor] = pre_residual[clusterNum][g[*i].int_descriptor] + tmp;
        //}
    }
}

void single_pagerank(graph& g, int clusterNum){
    vertex_iterator i,j;

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
    double RxSum, RySum;
    for(int q = 0; q < rankiter; q++){
        RxSum = 0;
        RySum = 0;
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
            if(g[*i].int_descriptor < xNum){
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