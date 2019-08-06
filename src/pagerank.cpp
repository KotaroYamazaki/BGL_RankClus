
#include "pagerank.hpp"
#include "gauss_southwell.hpp"

pagerank::pagerank(int _xNum, int _t, int _clusterNum, double _epsi){
    xNum = _xNum;
	t = _t;
    epsi = _epsi;
	clusterNum = _clusterNum;
}

void pagerank::solve(graph& g, int iteration_num){
    normalize_outedge_weight(g);
    if(iteration_num == 0){
        pagerank_from_scratch(g);
        get_rank_for_rankclus(g);
    }else{
        if(t < gauss_start){
            pagerank_from_scratch(g);
            get_rank_for_rankclus(g);
            gauss_southwell gs(epsi, alpha, clusterNum);
            gs.init(g);
        }else{
            gauss_southwell gs(epsi, alpha, clusterNum);
            gs.solve(g);
            get_rank_for_rankclus(g);
            gs.update_pregraph(g);
        }   
    }
}

void pagerank::normalize_outedge_weight(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
            for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                g[*e .first].weight /= out_degree(*i, g);
            }
    }
}

void pagerank::pagerank_from_scratch(graph& g){
    // bias value
    double b = 1.0/num_vertices(g);
    vector<double> rank;
    rank = vector<double>(num_vertices(g), b);
    vector<double> tmp_rank = rank;
    
    vertex_iterator i,j;
    int v;
    bool conv_flag = false;
    //for(v = 0; v < rankiter; v++){
    while(!conv_flag){
	conv_flag = true;
	double change = 0;
        double ranksum = 0;

        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                tmp_rank[*i] += g[*e.first].weight * rank[source(*e.first, g)];
            }
            tmp_rank[*i] = alpha * tmp_rank[*i] + (1 - alpha) * b;
            ranksum += tmp_rank[*i];
        }

        for(boost::tie(i,j) = vertices(g); i !=j; i++){
            tmp_rank[*i] /= ranksum;

            change = fabs(rank[*i] - tmp_rank[*i]);
            if(change > epsi)conv_flag = false;
            
            rank[*i] = tmp_rank[*i];
            g[*i].p_rank = rank[*i];
            tmp_rank[*i] = 0;
        }
        if(conv_flag)break;
    }
}

void pagerank::get_rank_for_rankclus(graph& g){
    //Calc rankscore for rankclus.(Convert single-graph rank score -> bi-type network graph.)
    vertex_iterator i, j;
    double rxsum = 0;
    double rysum = 0;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(clusterNum)){
            rxsum += g[*i].p_rank;
        }else{
            rysum += g[*i].p_rank;
        }
    }

    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        if(*i < xNum && g[*i].belongs_to_cluster(clusterNum)){
            g[*i].rx = g[*i].p_rank/rxsum;
        }else{
            g[*i].ry = g[*i].p_rank/rysum;
        }
    }
}
