#include"clustering.hpp"

void print_cluster_with_label(graph& g);
bool has_empty_cluster(graph& g);

clustering::clustering(vector<double>& _WkXY_sum, int _K, double _WXY_sum,int _xNum, vector<vector<int>>& _cluster_label){
    WkXY_sum = _WkXY_sum;
    K = _K;
    WXY_sum = _WXY_sum;
    xNum = _xNum;
    cluster_label = _cluster_label;
}

double clustering::Norm(vector<double>& array ){
	double Sum = 0;
	for(int l = 0;l < K; l++){
		Sum += array[l]*array[l];
	}
	return(sqrt(Sum));
}

vector<double> clustering::get_generating_probability(){
    vector<double> p;
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        p.push_back(1.0*WkXY_sum[clusterNum]/WXY_sum);
    }
    return p;
}
vector<double> clustering::calc_conditional_distribution(const graph& g, const vector<graph>& subgraph,const vector<double> p){
        vertex_iterator i,j;
        vector<double> conditional_p;
        for(int z = 0; z < K; z++){
            double tmp_p = 0;
            for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
                for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    //if(isnan(subgraph[z][source(*e.first, g)].ry))cout << subgraph[z][source(*e.first, g)].ry << endl;
                    tmp_p += g[*e.first].weight * subgraph[z][target(*e.first, g)].rx *subgraph[z][source(*e.first, g)].ry * p[z];
                }
        }
        conditional_p.push_back(tmp_p/WXY_sum);
        //if(isnan(p[z]))cout << "p[" << z << "]:"<<  p[z] << endl;
    }
    return conditional_p;
}

vector<vector<double>> clustering::calc_pi_using_bayesian_rule(const graph& g, const vector<graph>& subgraph,const vector<double> p){
    vertex_iterator i,j;
    vector<vector<double>> pi = vector<vector<double>>(K);
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
        }
    }
    return pi;
}

vector<vector<double>> clustering::get_K_dimentional_vector(const vector<vector<double>>& pi){
    vector<vector<double>> s = vector<vector<double>>(xNum);
    for(int i = 0; i < xNum; i++){
        for(int k = 0; k < K; k++){
            s[i].push_back(pi[k][i]);
        }
    }
    return s;
}

vector<vector<double>> clustering::get_center_vector(const graph& g, const vector<vector<double>>& s){
    vector<vector<double>> center_vec;
    center_vec = vector<vector<double>>(K,vector<double>(K,0));

    for( int Xk = 0; Xk < K; Xk++){
        int cluster_size = cluster_label[Xk].size();
        for(int col = 0; col < K; col++){
            vertex_iterator i,j;
            for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
                if(g[*i].cluster_label == Xk){
                    center_vec[Xk][col] += s[*i][col];
                }
            }
            center_vec[Xk][col] /= cluster_size;
        }
    }
    return center_vec;
}

double clustering::calc_distance_by_one_minus_cosine_similarity(vector<double> s_x, vector<vector<double>>center_vec, int k){
    double tmp = 0;
    for (int l = 0; l < K; l++){
        tmp += s_x[l] * center_vec[k][l];
    }
    double cosine_similarity = 1.0 - (tmp/(Norm(s_x) * Norm(center_vec[k])));
    return cosine_similarity;
}

int clustering::get_index_of_nearest_cluster(vector<double> s_x,const vector<vector<double>> center_vec){
    vector<double> D;
        int index = -1;
        double minDis = 10;
        for (int k = 0; k < K; k++){
            D.push_back(calc_distance_by_one_minus_cosine_similarity(s_x,center_vec, k));
            if(D[k] < minDis){
                minDis = D[k];
                index = k;
            }
        }
    return index;
}

void clustering::update_cluster_label(vertex_property& v, const int index){
    if(v.cluster_label == index){
            v.same_previous_cluster = true;
        }else{
            v.same_previous_cluster = false;
        }
        v.cluster_label = index;
    }

void clustering::check_empty_cluster(graph & g){
    if(has_empty_cluster(g)){
        cout << "Cluster became empty , you have to run the program again." << endl;
        print_cluster_with_label(g);
        exit(0);
    }
}

vector<vector<int>> clustering::update_cluster_label(graph &g, vector<graph>& subgraph){
    auto p = get_generating_probability();
    auto conditional_p = calc_conditional_distribution(g, subgraph, p);
    auto pi = calc_pi_using_bayesian_rule(g, subgraph, conditional_p);

    auto s = get_K_dimentional_vector(pi);
    auto center_vec = get_center_vector(g, s);
    vertex_iterator m,n;
    vector<vector<int>> new_cluster_label(K);
    for (boost::tie(m, n) = vertices(g); *m < xNum; m++) {
        int index = get_index_of_nearest_cluster(s[*m], center_vec);
        update_cluster_label(g[*m], index);
        new_cluster_label[index].push_back(*m);
    }
    return new_cluster_label;
}