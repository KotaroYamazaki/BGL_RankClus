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

double Norm(vector<double>& array ){
	double Sum = 0;
	for(int l = 0;l < K; l++){
		Sum += array[l]*array[l];
	}
	return(sqrt(Sum));
}

graph clustering(graph &g, vector<graph>& subgraph){
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
                //サブグラフでエッジないとこにあくせすしてしまう？　
                //それは違うか
                tmp_p += g[*e.first].weight * subgraph[z][target(*e.first, g)].rx *subgraph[z][source(*e.first, g)].ry * p[z];
            }
        }
        p[z] = tmp_p/WXY_sum;
    }

    // calc pi
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
        double tmp_sum = 0;
        for(int l = 0; l < K; l++){
            tmp_sum += subgraph[l][*i].rx * p[l];
        }
        //cout << "tmp_dum: " << tmp_sum << endl;

        for(int z = 0; z < K; z++){ 
            //cout << "condirank[" << z << "](" << *i << "): " << subgraph[z][*i].conditional_rank *p[z] << endl;
            double val = subgraph[z][*i].conditional_rank *p[z]/tmp_sum;
           // if(isnan(subgraph[z][*i].conditional_rank)) cout << "nandesu " << endl;
            pi[z].push_back(val);
            //cout << val << endl;
        }
    }

    vector<vector<double>> s(xNum);
    for(int i = 0; i < xNum; i++){
        for(int k = 0; k < K; k++){
            s[i].push_back(pi[k][i]);
        }
    }
    
    vector<vector<double>> center_vec(K);
    for (int dim = 0; dim < K; dim++){
        for(int clusterNum = 0; clusterNum < K; clusterNum++){
            double tmp_center = 0;
            for (boost::tie(i, j) = vertices(g); *i < xNum; i++) {
                if(g[*i].belongs_to_cluster == clusterNum){
                    tmp_center += s[*i][clusterNum];
                }
            }
            center_vec[dim].push_back(tmp_center);
        }
    }
    
    vector<vector<double> > D;
	//比較用前のラベル

	D = vector<vector<double>>(xNum,vector<double>(K,0));
	//ラベルの初期化
		for(int i = 0; i < xNum; i++){
		//Dの計算
			double normpi = Norm(pi[i]);
            if(normpi == 0) cout << "odayo" << endl;
			int index = 0;
			double minDis = 1;
			for (int k = 0; k < K; k++){
				double tmp = 0;
				for (int l = 0; l < K; l++){
					tmp += pi[i][l]*center_vec[k][l];
				}
                //cout << "check" << endl;
                double tmp2 = (normpi*Norm(center_vec[k]));
                if(tmp2 == 0) cout << "odayo" << endl;
                D[i][k] = tmp2;
		//assign
				// if(D[i][k] < minDis){
				// 	minDis = D[i][k];
				// 	index = k;
				// }
			}

		}
	
    
    return g;
}

