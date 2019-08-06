#include "gauss_southwell.hpp"

vector<graph> pre_graph;
vector<vector<double>> residual;

gauss_southwell::gauss_southwell(double _epsi, double _alpha,int _clusterNum){
	epsi = _epsi;
	alpha = _alpha;
	clusterNum = _clusterNum;
}

void gauss_southwell::solve(graph &g){
	pair<queue<int>, vector<bool>> p = calc_tracking_residual(g);
	queue<int> q = p.first;
	vector<bool> occupied_flag = p.second;
	vector<double> &res = residual[clusterNum];

	while (!q.empty())
	{
		unsigned long index = q.front();
		vertex_descriptor max_index = vertex(index, g);
		q.pop();
		occupied_flag[max_index] = false;

		double r_i = res[max_index];
		//x^(v) = x^(v-1) + r_i^(v-1)*e_i
		g[max_index].p_rank += r_i;

		//r^(v) = r^(v-1) - r_i^(v-1)*e_i
		res[max_index] -= r_i;

		//r^(v) = r^(v-1) + a*r_i^(v-1)*P*e_i
		for (auto e = out_edges(max_index, g); e.first != e.second; e.first++)
		{
			unsigned long index = target(*e.first, g);
			res[index] += alpha * (g[*e.first].weight * r_i);
			if (fabs(res[index]) > epsi && !occupied_flag[index])
			{
				q.push(index);
				occupied_flag[index] = true;
			}
		}
	}
}

void gauss_southwell::update_pregraph(graph& g){
	pre_graph[clusterNum] = g;
}

void gauss_southwell::init(graph& g){
	calc_initial_residual(g);
	init_pregraph(g);
}

void gauss_southwell::init_pregraph(graph& g){
	pre_graph.push_back(g);
}

void gauss_southwell::calc_initial_residual(graph &g)
{
	vector<double> tmp_res;
	double b = 1.0 / num_vertices(g);
	vertex_iterator i, j;
	for (boost::tie(i, j) = vertices(g); i != j; i++)
	{
		double tmp = (1 - alpha) * b;
		for (auto e = in_edges(*i, g); e.first != e.second; e.first++)
		{
			tmp += alpha * g[*e.first].weight * g[source(*e.first, g)].p_rank;
		}
		tmp -= g[*i].p_rank;
		tmp_res.push_back(tmp);
	}

	residual.push_back(tmp_res);


}

pair<queue<int>, vector<bool>> gauss_southwell::calc_tracking_residual(graph &g)
{
	queue<int> q;
	vector<bool> occupied_flag = vector<bool>(num_vertices(g));
	vertex_iterator i, j;
	graph &pre_g = pre_graph[clusterNum];
	for (boost::tie(i, j) = vertices(g); i != j; i++)
	{
		g[*i].p_rank = pre_g[*i].p_rank;
		vertex_descriptor v;
		double Pt_times_x = 0;
		double Pt_mius_1_times_x = 0;
		//calc P(t)*x(t-1)
		for (auto e = in_edges(*i, g); e.first != e.second; e.first++)
		{
			v = source(*e.first, g);
			Pt_times_x += g[*e.first].weight * pre_g[v].p_rank;
		}
		//calc P(t-1)*x(t-1)
		for (auto e = in_edges(*i, pre_g); e.first != e.second; e.first++)
		{
			v = source(*e.first, pre_g);
			Pt_mius_1_times_x += pre_g[*e.first].weight * pre_g[v].p_rank;
		}

		residual[clusterNum][*i] += alpha * (Pt_times_x - Pt_mius_1_times_x);
		if (fabs(residual[clusterNum][*i]) > epsi)
		{
			q.push(*i);
			occupied_flag[*i] = true;
		}
	}
	return make_pair(q, occupied_flag);
}
