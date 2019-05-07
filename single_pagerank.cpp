#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <time.h>
using namespace std;
extern vector<double> WkXY_sum;
extern int xNum;

const double alpha = 0.85;
const int rankiter = 15;
const int gauss_itr = 20000;
const double epsi = 0.01;

extern int K;
vector<graph> pre_graph;
vector<vector<double>> residual;
vector<vector<double>> pre_residual;
void gauss_southwell(graph& g, int clusterNum);
void normalize_outedge_weight(graph& g);
void normalize_xy_rank(graph& g, int clusterNum);
void normalize_global_rank(graph& g);

extern int t;
extern int iteration_num;
extern vector<vector<int>> cluster_label;

void single_pagerank(graph& g, int clusterNum);
void authority_ranking(graph& g, int clusterNum);
vector<double> init_rank(graph& g);
void init_residual(graph& g, int clusterNum);

void ranking(graph& subgraph, int clusterNum){
    //clock_t start = clock();
    if(iteration_num == 1){
        normalize_outedge_weight(subgraph);
        if(t == 0){
            //authority_ranking(subgraph, clusterNum);
            single_pagerank(subgraph, clusterNum);
            pre_graph.push_back(subgraph);
        }else{
            normalize_global_rank(subgraph);
            init_residual(subgraph, clusterNum);
            gauss_southwell(subgraph, clusterNum);
            normalize_xy_rank(subgraph, clusterNum);
            pre_graph[clusterNum] = subgraph;
        }
    }else{
        authority_ranking(subgraph, clusterNum);
    }
}



void normalize_outedge_weight(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        //O(n*E*E)
        //W_XY で行正規化をしている
        //if(g[*i].int_descriptor > xNum){
            //g[*i].rx = pre_graph[clusterNum][*i].rx;
            int rowsum = 0;
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                rowsum += g[*e .first].weight;
            }
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                g[*e .first].weight /= rowsum;
                //cout << "out wieight : " << g[*e.first].weight << endl;
            }
     //   }
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
                g[*i].rx /= RxSum;
                // cout<< g[*i].rx << endl;
            }else{
                g[*i].ry /= RySum;
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


vertex_iterator cast_vertex_iterator(unsigned long index, graph g){
    vertex_iterator i,j;
    for(boost::tie(i,j) = vertices(g); g[*i].int_descriptor < index; i++){
    }
    return i;
}

void gauss_southwell(graph& g, int clusterNum){
    for(int v = 0; v < gauss_itr; v++){
        auto max_itr = max_element(residual[clusterNum].begin(), residual[clusterNum].end());
        unsigned long max_index = distance(residual[clusterNum].begin(), max_itr);
        double r_i = residual[clusterNum][max_index];
        vertex_iterator max_vertex_itr = cast_vertex_iterator(max_index, g);
        //if()file << r_i << ", ";
        // if(v > 0){
             //cout << "[" << v << "] r_i: " << r_i << " | index " << max_index << endl;
        //     //cout <<  "max_ri"<< residual[clusterNum][max_index] << endl;
        //     //cast_vertex_iterator(max_index, g);
        //     // cout << *max_vertex_itr << endl;
        //     // cout << "rx" << g[*max_vertex_itr].ry << endl;
        // }

        if(max_index < xNum){
            g[*max_vertex_itr].rx += r_i;
        }else{
            g[*max_vertex_itr].ry += r_i;
        }
        residual[clusterNum][max_index] -= r_i;

        for (auto e = out_edges(*max_vertex_itr, g); e.first!=e.second; e.first++) {
            //cout << residual[clusterNum][target(*e.first, g)] << endl;
           residual[clusterNum][target(*e.first, g)] +=  alpha * (g[*e.first].weight * r_i);
            //cout << "伝搬後　"<<  alpha << "*"<< (g[*e.first].weight) << "*" <<  r_i << " = " << residual[clusterNum][target(*e.first, g)] << endl;
        }

        if(r_i < epsi || v == gauss_itr - 1){
            cout << "converge at " << v + 1<< endl;
            break;
        }
    }
}

void init_residual(graph& g, int clusterNum){
    if(t >= 1){
        residual = vector<vector<double>>(K,vector<double>(num_vertices(g),0));
    }
    

        vertex_iterator i,j, m, n;
        for(boost::tie(i,j) = vertices(g); i != j; i++){

            // ランク値の引き継ぎ
            if(g[*i].int_descriptor < xNum){
                g[*i].rx = pre_graph[clusterNum][*i].rx;
                //cout << pre_graph[clusterNum][*m].rx << endl;
            }else{
                g[*i].ry = pre_graph[clusterNum][*i].ry;
            }

            double tmp_res = 0;
            for(boost::tie(m,n) = vertices(g); m != n; m++){
                double rank = 0;
                
                if(g[*m].int_descriptor < xNum)rank = pre_graph[clusterNum][*m].rx;
                else rank = pre_graph[clusterNum][*m].ry;
                // 前エッジあり　＆＆　今もある
                if(boost::edge(*i,*m, pre_graph[clusterNum]).second && boost::edge(*i, *m, g).second){
                    tmp_res += alpha * ( pre_graph[clusterNum][boost::edge(*i,*m, pre_graph[clusterNum]).first].weight - g[boost::edge(*i,*m, g).first].weight ) * rank;
                }else if(boost::edge(*i,*m, pre_graph[clusterNum]).second && !boost::edge(*i, *m, g).second){
                    //前あり　＆＆　今はない
                    tmp_res = alpha * ( pre_graph[clusterNum][boost::edge(*i,*m, pre_graph[clusterNum]).first].weight - 0 ) * rank;
                }else if(!boost::edge(*i,*m, pre_graph[clusterNum]).second && boost::edge(*i, *m, g).second){
                    //　まえない　＆＆　今ある
                    tmp_res = alpha * ( 0 - g[boost::edge(*i,*m, g).first].weight ) * rank;
                }else{
                    // ない　＆＆　ない
                    continue;
                }
            }
            residual[clusterNum][*i] += tmp_res;
        }
}

void single_pagerank(graph& g, int clusterNum){
    vertex_iterator i,j;
    double RxSum = 0;
    double RySum = 0;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(g[*i].int_descriptor < xNum)g[*i].rx = 1;
        else g[*i].ry = 1;

        for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
            if(g[*i].int_descriptor < xNum){
                g[*i].rx =  alpha * g[*e.first].weight * g[*i].rx + (1 - alpha)*(1/(cluster_label[clusterNum].size()));
                RxSum += g[*i].rx;
            }else{
                g[*i].ry = alpha * g[*e.first].weight;
                RySum += g[*i].ry;
            }
        }

    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            if(g[*i].int_descriptor < xNum){
                g[*i].rx /= RxSum;
                //cout << g[*i].rx << endl;
            }else{
                g[*i].ry /= RySum;
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