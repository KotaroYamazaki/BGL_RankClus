#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;
extern vector<double> WkXY_sum;
extern int K;
extern double WXY_sum;
extern int xNum;
bool has_empty_cluster(graph& g);
extern vector<vector<int>> cluster_label;

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
        //if(isnan(p[clusterNum]))cout << "p[" << clusterNum << "]:"<<  p[clusterNum] << endl;
    }

    for(int z = 0; z < K; z++){
        double tmp_p = 0;
        for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                //if(isnan(subgraph[z][source(*e.first, g)].ry))cout << subgraph[z][source(*e.first, g)].ry << endl;
                tmp_p += g[*e.first].weight * subgraph[z][target(*e.first, g)].rx *subgraph[z][source(*e.first, g)].ry * p[z];
            }
        }
        p[z] = tmp_p/WXY_sum;
        //if(isnan(p[z]))cout << "p[" << z << "]:"<<  p[z] << endl;
    }

    // calc pi
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        double tmp_sum = 0;
        for(int l = 0; l < K; l++){
            tmp_sum += subgraph[l][*i].conditional_rank * p[l];
            if(isnan(subgraph[l][*i].conditional_rank)){
                cout << "subgraph[l][*i].conditional_rank is nan";
                exit(1);
            }
        }
        for(int z = 0; z < K; z++){ 
            double val = 0;
            if(tmp_sum != 0)val = subgraph[z][*i].conditional_rank * p[z]/tmp_sum;
            pi[z].push_back(val);
            //cout << "pi[" << z << "]:"<<  val << endl;
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

    for( int Xk = 0; Xk < K; Xk++){
        int cluster_size = cluster_label[Xk].size();
        for(int col = 0; col < K; col++){
            for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
                if(g[*i].cluster_label == Xk){
                    center_vec[Xk][col] += s[*i][col];
                }
            }
            center_vec[Xk][col] /= cluster_size;
        }
    }
    
    // 距離を格納するリスト
    vector<vector<double> > D;
	D = vector<vector<double>>(xNum,vector<double>(K,0));
    vertex_iterator m,n;

    vector<vector<int>> new_cluster_label(K);
        
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
            //cout <<  *m << "->" << k << ": " << D[*m][k] << endl;
            if(D[*m][k] < minDis){
                minDis = D[*m][k];
                index = k;
            }
        }

        if(g[*m].cluster_label == index){
            g[*m].same_previous_cluster = true;
        }else{
            g[*m].same_previous_cluster = false;
        }
        g[*m].cluster_label = index;
        new_cluster_label[index].push_back(*m);
    }
    cluster_label = new_cluster_label;

    // for(int i = 0; i < K; i++){
    //     cout << "cluster = { ";
    //     for(int j = 0; j < cluster_label[i].size(); j++){
    //         cout << cluster_label[i][j] << ", "<< flush;
    //     }
    //     cout << "}" <<endl;
    // }

    void print_cluster_with_label(graph& g);
    if(has_empty_cluster(g)){
        cout << "Cluster became empty , you have to run the program again." << endl;
        print_cluster_with_label(g);
        exit(0);
    }
}