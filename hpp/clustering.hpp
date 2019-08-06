#ifndef _CLUSTERING_H_
#define _CLUSTERING_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

class clustering{
    private:
        vector<double> WkXY_sum;
        int K;
        double WXY_sum;
        int xNum;
        vector<vector<int>> cluster_label;
        
        double Norm(vector<double>& array );
        vector<double> get_generating_probability();
        vector<double> calc_conditional_distribution(const graph& g, const vector<graph>& subgraph,const vector<double> p);
        vector<vector<double>> calc_pi_using_bayesian_rule(const graph& g, const vector<graph>& subgraph,const vector<double> p);
        vector<vector<double>> get_K_dimentional_vector(const vector<vector<double>>& pi);
        vector<vector<double>> get_center_vector(const graph& g, const vector<vector<double>>& s);
        double calc_distance_by_one_minus_cosine_similarity(vector<double> s_x, vector<vector<double>>center_vec, int k);
        int get_index_of_nearest_cluster(vector<double> s,const vector<vector<double>> center_vec);
        void update_cluster_label(vertex_property& v, const int index);
        void check_empty_cluster(graph & g);

    public:
        clustering(vector<double>& _WkXY_sum, int _K, double _WXY_sum,int _xNum, vector<vector<int>>& _cluster_label);
        vector<vector<int>> update_cluster_label(graph &g, vector<graph>& subgraph);
};

#endif