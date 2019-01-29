#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;
extern vector<int> WkXY_sum;
extern int K;
extern int WXY_sum;
extern int xNum;

vector<int> calc_cluster_size(graph& g){
    vector<int> size(K,0);
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
        size[g[*i].belongs_to_cluster] += 1;
    }
    return size;
}

double Norm(vector<double>& array ){
	double Sum = 0;
	for(int l = 0;l < K; l++){
		Sum += array[l]*array[l];
	}
	return(sqrt(Sum));
}

void clustering(graph &g, vector<graph>& subgraph){
    vector<vector<double>> pi(K);
    vector<double> p;
    vertex_iterator i,j;
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        p.push_back(1.0*WkXY_sum[clusterNum]/WXY_sum);
    }

    for(int z = 0; z < K; z++){
        
        double tmp_p = 0;
        for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                tmp_p += g[*e.first].weight * subgraph[z][target(*e.first, g)].rx *subgraph[z][source(*e.first, g)].ry * p[z];
            }
        }
        p[z] = tmp_p/WXY_sum;
    }

    // calc pi
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
        double tmp_sum = 0;
        for(int l = 0; l < K; l++){
            
            tmp_sum += subgraph[l][*i].conditional_rank * p[l];
        }
        for(int z = 0; z < K; z++){ 
            double val = 0;
            if(tmp_sum != 0)val = subgraph[z][*i].conditional_rank * p[z]/tmp_sum;
            pi[z].push_back(val);
        }
    }

    vector<vector<double>> s(xNum);
    for(int i = 0; i < xNum; i++){
        for(int k = 0; k < K; k++){
            s[i].push_back(pi[k][i]);
        }
    }
    
    vector<vector<double>> center_vec;
    center_vec = vector<vector<double>>(K,vector<double>(K,0));
    vector<int> cluster_size = calc_cluster_size(g);

    for( int Xk = 0; Xk < K; Xk++){
        for(int col = 0; col < K; col++){
            for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
                if(g[*i].belongs_to_cluster == Xk){
                    center_vec[Xk][col] += s[*i][col];
                }
            }
            center_vec[Xk][col] /= cluster_size[Xk];
        }
    }
    
    // 距離を格納するリスト
    vector<vector<double> > D;
	D = vector<vector<double>>(xNum,vector<double>(K,0));
        vertex_iterator m,n;

        for (boost::tie(m, n) = vertices(g); *m < xNum; m++) {
		//Dの計算
			double norms = Norm(s[*m]);
			int index = 0;
			double minDis = 1;
			for (int k = 0; k < K; k++){
				double tmp = 0;
				for (int l = 0; l < K; l++){
					tmp += s[*m][l] * center_vec[k][l];
				}
                D[*m][k] = 1.0 - (tmp/(norms * Norm(center_vec[k])));
		//assign
				if(D[*m][k] < minDis){
					minDis = D[*m][k];
					index = k;
				}
			}

            if(g[*m].belongs_to_cluster == index){
                g[*m].same_previous_cluster = true;
            }else{
                g[*m].same_previous_cluster = false;
            }
            g[*m].belongs_to_cluster = index;
		}
}

