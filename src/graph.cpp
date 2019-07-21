#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
#include <algorithm>
using namespace std;

int xNum;
int yNum;
double WXY_sum = 0;
vector<double> WkXY_sum;
extern int K;
int top_k = 10;
extern int input_seed;
extern string path;
extern vector<vector<int>> cluster_label;

vector<string> name_vector;
vector<vector<double>> row_sum_vec;

vector<string> split(string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

graph construct_graph(){
    // エッジのリスト
    vector<edge> edge_vector;
    // 各エッジの属性値の構造体のリスト
    vector<edge_property> property_vector;

    string str, name_X, name_Y;
    string file_WXY = path + "WXY.csv";
    string file_WYY = path + "WYY.csv";
    string file_X = path + "X.txt";
    string file_Y = path + "Y.txt";

    ifstream ifs_WXY(file_WXY);
    ifstream ifs_X(file_X);
    ifstream ifs_WYY(file_WYY);
    ifstream ifs_Y(file_Y);

    if(ifs_X.fail()){
            cout << "Error! Failed to read " << file_X  << "."<< endl;
            exit(0);
    }
    if(ifs_Y.fail()){
            cout << "Error! Failed to read " << file_Y  << "."<< endl;
            exit(0);
    }
    if(ifs_WXY.fail()){
            cout << "Error! Failed to read " << file_WXY  << "."<< endl;
            exit(0);
    }
    if(ifs_WYY.fail()){
            cout << "Warning: Execute RankClus without WYY file:" << file_WYY  << "."<< endl;
    }

    

    xNum = 0;
    while(getline(ifs_X, name_X)){
        name_vector.push_back(name_X);
        xNum += 1;
    }
    yNum = 0;
    while(getline(ifs_Y, name_Y)){
        name_vector.push_back(name_Y);
        yNum += 1;
    }

    WXY_sum = 0;
    int from, to, val;
    // Target type
	while(getline(ifs_WXY, str)){
        vector<string> strvec = split(str, ',');
        from = stoi(strvec.at(0));
        to = stoi(strvec.at(1));
        val = stoi(strvec.at(2));

        struct edge_property a;
        a.label = "AtoT";
        a.weight = val;
        property_vector.push_back(a);
        from += xNum;
        WXY_sum += val;
		edge_vector.push_back(edge(from, to));
        //逆方向
        a.label = "TtoA";
        property_vector.push_back(a);
        edge_vector.push_back(edge(to, from));
	}

    // Attribute Type
    while(getline(ifs_WYY,str)){
        vector<string> strvec = split(str, ',');
        from = stoi(strvec[0]);
        to = stoi(strvec[1]);
        val = stoi(strvec[2]);

        if(from != to){
            struct edge_property a;
            a.label = "AtoA";
            a.weight = val;
            property_vector.push_back(a);
            from += xNum;
            to += xNum;
            edge_vector.push_back(edge(from, to));
            //　逆方向
            a.label = "AtoA";
            property_vector.push_back(a);
            edge_vector.push_back(edge(to, from));
        }
	}

    // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
    auto tag = boost::edges_are_unsorted_multi_pass;
    // グラフのコンストラクタ
    // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
    graph g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), xNum + yNum);

    return g;
}

void get_intial_partitions(graph& g){
    //random_device rnd;
    mt19937 rnd(input_seed);
    vertex_iterator i, j;
    vector<bool> check_flag(K,false);
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++){
		int num = rnd() % K ;
		//if(label[tmp].size() < m/K)//偏らないように調整
		//{label[tmp].push_back(i);}
		//else {i -= 1;}
        g[*i].cluster_label = num;
        cluster_label[num].push_back(*i);
        check_flag[num] = true;
    }
    //check
    for(auto itr = check_flag.begin(); itr != check_flag.end(); itr++){
        if(check_flag[*itr] == false){
            get_intial_partitions(g);
        }
    }
}

bool has_empty_cluster(graph& g){
    vector<bool> check_flag(K,false);
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++){
        check_flag[g[*i].cluster_label] = true;
    }
    for(int i = 0; i < K; i++)if(check_flag[i] == false) return true;
    return false;
}

string cast_state(const vertex_trajectory& State){
    if(State == NOTHING)return "NOTHING";
    if(State == STAY)return "STAY";
    if(State == ADD)return "ADD";
    if(State == LEAVE)return "LEAVE";
    return "NULL";
}

void update_added_vertex_state_for_cluster(vertex_trajectory& State){
    if(State == NOTHING || State == LEAVE) State = ADD;
    else State = STAY;
}

void update_not_belongs_vertex_state_for_cluster(vertex_trajectory& State){
    if(State == NOTHING || State == LEAVE) State = NOTHING;
    else State = LEAVE;
}

int count_number_of_change_edges(const vertex_trajectory s, const vertex_descriptor v, const graph& g){
    if(s == LEAVE || s == ADD)return (out_degree(v,g) + in_degree(v, g));
    else return 0;
}

vector<graph> construct_sub_graph(graph& g){
    //　サブグラフを格納するリスト
    vector<graph> subgraph_vector;
    vector<string> x_name;
    
    WkXY_sum = vector<double>(K,0);
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        vertex_iterator i,j;
        int cluster_size = 0;
        double row_sum = 0;
        int max_degree = 0;
        // エッジのリスト
        std::vector<edge> edge_vector;
        // 各エッジの属性値の構造体のリスト
        std::vector<edge_property> property_vector;
        int count_edge_change = 0;
        int edge_sum = 0;
        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
                // ノードがターゲットタイプかつ該当クラスタに所属する場合以下の処理を行う
                    // アトリビュートタイプノードからはいってくるエッジのプロパティをコピー
                    // アトリビュートタイプへでていくエッジのプロパティをコピー
            if(g[*i].int_descriptor < xNum){ 
                if(g[*i].belongs_to_cluster(clusterNum)){
                    cluster_size += 1;
                    x_name.push_back(g[*i].name);
                    for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                    // target(edge_descriptor, graph) でエッジの先のノードを取得できる
                    // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
                    // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
                        struct edge_property a;
                        a.label = "TtoA";
                        a.weight = g[*e.first].weight;
                        property_vector.push_back(a);
                        edge_vector.push_back(edge(source(*e.first, g), target(*e.first, g)));
                        edge_sum += a.weight;
                    }
                    update_added_vertex_state_for_cluster(g[*i].state[clusterNum]);
                }else{
                    update_not_belongs_vertex_state_for_cluster(g[*i].state[clusterNum]);
                }
                if(g[*i].state[clusterNum] == ADD || g[*i].state[clusterNum] == LEAVE)max_degree = max(max_degree, (int)(in_degree(*i, g)+out_degree(*i, g)));
                count_edge_change += count_number_of_change_edges(g[*i].state[clusterNum], *i, g);
                //cout << out_degree(*i , g) + in_degree(*i, g)<< endl;
                //if(count_edge_change != 0)cout << "count : "<< count_edge_change << "[" << max_degree << "]"<< endl;
            } else {
                for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
                    //　入次してくるエッジのもとのノードのタイプがattributeのときに以下の処理を行う（ターゲットタイプには行わない）
                    struct edge_property a;
                    if(g[target(*e.first,  g)].label == "attribute"){
                            a.label = "AtoA";
                            a.weight = g[*e.first].weight;
                            property_vector.push_back(a);
                            edge_vector.push_back(edge(source(*e.first, g), target(*e.first, g)));
                            // 逆方向
                            property_vector.push_back(a);
                            edge_vector.push_back(edge(target(*e.first, g), source(*e.first, g)));
                    }else{
                        if(g[target(*e.first, g)].cluster_label == clusterNum){
                            a.label = "AtoT";
                            a.weight = g[*e.first].weight;
                            property_vector.push_back(a);
                            edge_vector.push_back(edge(source(*e.first, g), target(*e.first, g)));
                        }
                    }
                }
            }
        }
        WkXY_sum[clusterNum] = edge_sum;
        // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
        auto tag = boost::edges_are_unsorted_multi_pass;
        // グラフのコンストラクタ
        // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
        graph sub_g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), xNum + yNum);
        subgraph_vector.push_back(sub_g);
        //cout << "global :: "<< num_edges(sub_g) << endl;
        }
    return subgraph_vector;
}

// グラフの初期化
void init_graph(graph& g){
    // ノードの走査
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        // (*vertex_iterator) で vertex_descriptor になる
        // g[vertex_descriptor].属性（vertex_property で定義したもの）で参照できる
        //name の入れ方name_vectorがグローバル変数
        if(g[*i].int_descriptor < xNum ){
            g[*i].label = "target";
            g[*i].rx = 0;
            g[*i].name = name_vector[*i];
            g[*i].conditional_rank = 0;
            g[*i].cluster_label = -1;
            g[*i].state = vector<vertex_trajectory>(K, NOTHING);
        } else {
            g[*i].label = "attribute";
            g[*i].ry = 0;
            g[*i].name = name_vector[*i];
        }
        g[*i].int_descriptor = static_cast<int>(*i);
    }
}

void init_graph(graph& g, graph& global_g){
    // ノードの走査
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        // (*vertex_iterator) で vertex_descriptor になる
        // g[vertex_descriptor].属性（vertex_property で定義したもの）で参照できる
        //name の入れ方name_vectorがグローバル変数
        if(g[*i].int_descriptor < xNum ){
            g[*i].label = "target";
            g[*i].rx = 0;
            g[*i].name = name_vector[*i];
            g[*i].cluster_label = global_g[*i].cluster_label;
        } else {
            g[*i].label = "attribute";
            g[*i].ry = 0;
            g[*i].name = name_vector[*i];
        }
        g[*i].int_descriptor = static_cast<int>(*i);
    }
}

void print_graph_detail(graph& g){
    cout << "# vertices: " << num_vertices(g) << endl;
    cout << "# edges   : " << num_edges(g) << endl;

    vertex_iterator i, j;
    cout << " --- adjacent vertices --- " << endl;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        cout <<  *i << " --> " << flush;
        // adjacent_vertices(vertex_descriptor, graph) で隣接ノードを取得できる
        for (auto itr = adjacent_vertices(*i, g); itr.first!=itr.second; itr.first++) {
            cout << *itr.first << ", " << flush;
        }
        cout << endl;
    }


    // out_edges(vertex_descriptor, graph) でノードの出エッジが取得できる
    cout << " --- edges --- " << endl;
    for (boost::tie(i, j) = vertices(g); i!=j; i++) {
        for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
            // target(edge_descriptor, graph) でエッジの先のノードを取得できる
            // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
            // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
            cout << *i << " --> " << target(*e.first, g) << ": " << g[*e.first].weight << " (" << g[*e.first].label << ")" << endl;
        }
    }
}

vector<int> get_sorted_list(graph& g, int clusterNum){
    struct name_rank {
		double rank;
		int id;
		// 最後のconstを忘れると"instantiated from here"というエラーが出てコンパイルできないので注意
		bool operator<( const name_rank& right ) const {
			return rank < right.rank ;
		}
	};
    vector<name_rank> ranking_list;

    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor <xNum ; i++) {
        name_rank a;
        if(g[*i].belongs_to_cluster(clusterNum)){
            a.rank = g[*i].rx;
            a.id = *i;
            ranking_list.push_back(a);
        }
    }
    sort(ranking_list.rbegin(), ranking_list.rend());

    vector<int> sorted_id_list;
    for(int i = 0; i < ranking_list.size(); i++){
        sorted_id_list.push_back(ranking_list[i].id);
    }
    return sorted_id_list;
}
void print_rank_within_cluster(graph& g, int clusterNum){

    vector<int> sorted_id_list = get_sorted_list(g, clusterNum);
    cout << "<Cluster = " << clusterNum + 1<< "> (size: " << sorted_id_list.size() << ")" << endl;
    //top_kがクラスタのサイズより大きかったら変更
    int top = top_k;
    if(top_k > sorted_id_list.size())top = sorted_id_list.size();
	
    for(int i = 0; i < top; i++){
        cout << i +1 << ": " << g[sorted_id_list[i]].name << " ... ["<< g[sorted_id_list[i]].rx  << "]" <<endl;
	}
	cout << endl;
}

void print_cluster_with_name(graph& g ){
    vector<vector<string>> cluster(K);
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor <xNum ; i++) {
        cluster[g[*i].cluster_label].push_back(g[*i].name);
    }
	for (int k = 0; k < K; k++){
		cout << "  Cluster[" << k+1 << "] = { " << flush;
		for (int i = 0; i < cluster[k].size(); i++){
			cout << cluster[k][i] << flush;
			if(i != cluster[k].size()-1)cout << ", "<< flush;
		}
		cout << " }" << endl;
	}
}

void print_cluster_with_label(graph& g ){
	for (int k = 0; k < K; k++){
		cout << "  Cluster[" << k+1 << "] = { " << flush;
		for (int i = 0; i < cluster_label[k].size(); i++){
			cout << cluster_label[k][i] << "["<< cast_state(g[vertex(cluster_label[k][i], g)].state[k]) << "]"<<flush;
			if(i != cluster_label[k].size()-1)cout << ", "<< flush;
		}
		cout << " }" << endl;
	}
}


void write_result(vector<graph>& sub_g, string out_file){
    fstream file1, file2, file3;
    file1.open("result_clutsering_" + out_file, ios::out);
    file2.open("result_ranking_" + out_file, ios::out);
    file3.open("result_clustering_id_" + out_file, ios::out);

    for (int k = 0; k < K; k++){
    vector<int> sorted_id_list = get_sorted_list(sub_g[k], k);
		file1 << "###### Cluster" << k+1 << " (" << sorted_id_list.size() << ") ######" << endl;
        file2 << "###### Cluster" << k+1 << " (" << sorted_id_list.size() << ") ######" << endl;
        file3 << "cluster_label[" << k << "] = {" << flush; 
		for (int i = 0; i < sorted_id_list.size(); i++){
			file1 << sub_g[k][sorted_id_list[i]].name << endl;
            file2 << i + 1 << ": " << sub_g[k][sorted_id_list[i]].name << " ... ["<< sub_g[k][sorted_id_list[i]].rx  << "]" <<endl;
            file3 << sorted_id_list[i] << flush;
            if(i != sorted_id_list.size()-1)file3 << ", "<< flush;
		}
        file1 << endl;
        file2 << endl;
        file3 << " }" << endl;
	}
    cout << "==== Written the result of clustering to " << out_file << " === " <<  endl;
    file1.close();
    file2.close();
    file3.close();
}

void conditional_ranking(graph& g, graph& subgraph){
    vertex_iterator i,j;
    double ranksum = 0;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        double tmp = 0;
        for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
            tmp += g[*e.first].weight/out_degree(source(*e.first, g), g) * subgraph[source(*e.first, g)].ry; 
        }
        if(isnan(tmp)){
            cout << "tpm is nan" << endl;
            exit(1);
        }
        subgraph[*i].conditional_rank = tmp;
        ranksum += tmp;
    }
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++) {
        if(ranksum != 0)subgraph[*i].conditional_rank /= ranksum;
        else subgraph[*i].conditional_rank = 0;
        //cout << subgraph[*i].conditional_rank  << endl;
    }
}

bool check_converge_cluster(graph& g){
    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum ; i++) {
        if(!g[*i].same_previous_cluster) return false;
    }
    cout << "---------- CLUSTERS HAVE BEEN CONVERGED ---------" << endl;
    return true;
}

void write_result_to_csv(vector<int> time){
        ofstream file1;
        file1.open("results/result_time_compare.csv",ios_base::app);
        int comp_num = 2;	
            for(int i = 0; i< comp_num;i++){
            file1 << time[i] << flush;
            if(i != comp_num - 1)file1  << "," << flush;
            }
        file1.close();
    }

void write_result_for_NMI(graph& g, int iteration_num){
        string filename;
        if(iteration_num == 0)filename = "results/correct.csv";
        if(iteration_num == 1)filename = "results/result.csv";

        fstream file;
        file.open(filename,ios::out);

        vertex_iterator i,j;
        for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum ; i++) {
            file << g[*i].cluster_label << flush;
            if(g[*i].int_descriptor < xNum - 1) file << "," << flush;
        }
    }

    void write_result_in_each_iteration(graph& g, int iteration_num, int t){
        string filename;
        if(iteration_num == 0)filename = "results/cluster_result_rankclus.csv";
        if(iteration_num == 1)filename = "results/cluster_result_proposal.csv";

        fstream file;
        if(t == 0) file.open(filename, ios::out);
        else file.open(filename, ios::app);

        vertex_iterator i, j;
        file << "iteration" << t+1 << "," << flush;
        for (boost::tie(i, j) = vertices(g); g[*i].int_descriptor < xNum; i++)
        {
            file << g[*i].cluster_label << flush;
            if (g[*i].int_descriptor < xNum - 1)
                file << "," << flush;
        }
        file << endl;
        file.close();
    }