#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "graph.hpp"
using namespace std;

int xNum;
int yNum;
int WXY_sum =0;
vector<int> WkXY_sum;
extern int K;
int top_k = 10;

vector<string> name_vector;
vector<vector<string>> X_sub_name_vector;

vector<string> split(string& input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

ifstream read_file(string file_name){
    ifstream ifs(file_name);
        if(ifs.fail()){
            if(file_name != "fnameWYY"){
            cout << "Failed to read " << file_name  << "."<< endl;
            exit(0);
            }else{
                //WYY がなくても実行できるが警告を出す。
                cout << "Warnig: Excute with WYY file." << endl;
            }
        }
    return ifs;
}

graph construct_graph(){
    // エッジのリスト
    vector<edge> edge_vector;
    // 各エッジの属性値の構造体のリスト
    vector<edge_property> property_vector;

    string str, name_X, name_Y;
    string path = "dataset/dblp-small/";
    string fnameWXY = path + "WXY.csv";
    string fnameWYY = path + "WYY.csv";
    string file_X = path + "X.txt";
    string file_Y = path + "Y.txt";
    // string fnameWXY = "dataset/ml-latest/rating.csv";
    // string file_X = "dataset/ml-latest/movie.txt";
    // string file_Y = "dataset/ml-latest/user.txt";
    // string fnameWYY = "";

	ifstream ifs_WXY= read_file(fnameWXY);
    ifstream ifs_X = read_file(file_X);
    ifstream ifs_WYY = read_file(fnameWYY);
    ifstream ifs_Y = read_file(file_Y);

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

    //graph_property graph_prop;
    //graph_prop.edge_sum = Wsum;
    // tag は特に指定がなければ edges_are_unsorted_multi_pass で良い
    auto tag = boost::edges_are_unsorted_multi_pass;
    // グラフのコンストラクタ
    // エッジのコンテナの begin と end、エッジのプロパティのコンテナの begin、ノード数を渡す
    graph g(tag, edge_vector.begin(), edge_vector.end(), property_vector.begin(), xNum + yNum);

    return g;
}

void get_intial_partitions(graph& g){
    random_device rnd;
    vertex_iterator i, j;
    vector<bool> check_flag(K,false);
    for (boost::tie(i, j) = vertices(g); *i < xNum; i++){
		int num = rnd()%K;
		//if(label[tmp].size() < m/K)//偏らないように調整
		//{label[tmp].push_back(i);}
		//else {i -= 1;}
        g[*i].belongs_to_cluster = num;
        check_flag[num] = true;
    }

    //check
    for(auto itr = check_flag.begin(); itr != check_flag.end(); itr++){
        if(check_flag[*itr] == false){
            get_intial_partitions(g);
        }
    }
}

vector<graph> construct_sub_graph(graph& g){
    //　サブグラフを格納するリスト
    vector<graph> subgraph_vector;
    vector<string> x_name;
    X_sub_name_vector = vector<vector<string>>(K);
    
    WkXY_sum = vector<int>(K,0);
    for(int clusterNum = 0; clusterNum < K; clusterNum++){
        vertex_iterator i,j;
        int cluster_size = 0;
        // エッジのリスト
        std::vector<edge> edge_vector;
        // 各エッジの属性値の構造体のリスト
        std::vector<edge_property> property_vector;

        int edge_sum = 0;

        for (boost::tie(i, j) = vertices(g); i!=j; i++) {
                // ノードがターゲットタイプかつ該当クラスタに所属する場合以下の処理を行う
                    // アトリビュートタイプノードからはいってくるエッジのプロパティをコピー
                    // アトリビュートタイプへでていくエッジのプロパティをコピー
             if(g[*i].label == "target" && g[*i].belongs_to_cluster == clusterNum){ 
                    cluster_size += 1;
                    x_name.push_back(g[*i].name);
                    for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    // target(edge_descriptor, graph) でエッジの先のノードを取得できる
                    // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
                    // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
                        struct edge_property a;
                        a.label = "AtoT";
                        a.weight = g[*e.first].weight;
                        property_vector.push_back(a);
                        edge_vector.push_back(edge(source(*e.first, g), target(*e.first, g)));
                        //　逆方向
                        a.label = "TtoA";
                        property_vector.push_back(a);
                        edge_vector.push_back(edge(target(*e.first, g), source(*e.first, g)));
                        //エッジの重み合計
                        edge_sum += a.weight;
                    }
            }else if (g[*i].label == "attribute"){
                for (auto e = in_edges(*i, g); e.first!=e.second; e.first++) {
                    //　入次してくるエッジのもとのノードのタイプがattributeのときに以下の処理を行う（ターゲットタイプには行わない）
                    if(g[source(*e.first,  g)].label == "attribute"){
                            struct edge_property a;
                            a.label = "AtoA";
                            a.weight = g[*e.first].weight;
                            property_vector.push_back(a);
                            edge_vector.push_back(edge(source(*e.first, g), target(*e.first, g)));
                            // 逆方向
                            property_vector.push_back(a);
                            edge_vector.push_back(edge(target(*e.first, g), source(*e.first, g)));
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
        X_sub_name_vector.push_back(x_name);
        subgraph_vector.push_back(sub_g);
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
        if(*i < xNum ){
            g[*i].label = "target";
            g[*i].rx = 0;
            g[*i].name = name_vector[*i];
            g[*i].conditional_rank = 0;
            g[*i].belongs_to_cluster = -1;
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

    // vertex_iterator i, j;
    // cout << " --- adjacent vertices --- " << endl;
    // for (boost::tie(i, j) = vertices(g); i!=j; i++) {
    //     cout <<  *i << " --> " << flush;
    //     // adjacent_vertices(vertex_descriptor, graph) で隣接ノードを取得できる
    //     for (auto itr = adjacent_vertices(*i, g); itr.first!=itr.second; itr.first++) {
    //         cout << *itr.first << ", " << flush;
    //     }
    //     cout << endl;
    // }


    // // out_edges(vertex_descriptor, graph) でノードの出エッジが取得できる
    // cout << " --- edges --- " << endl;
    // for (boost::tie(i, j) = vertices(g); i!=j; i++) {
    //     for (auto e = out_edges(*i, g); e.first!=e.second; e.first++) {
    //         // target(edge_descriptor, graph) でエッジの先のノードを取得できる
    //         // out_edges() ではなく、in_edges() の場合は source() でエッジの元のノードを取得できる
    //         // g[edge_descriptor].属性 で、edge_property で定義した属性を参照できる
    //         cout << *i << " --> " << target(*e.first, g) << ": " << g[*e.first].weight << " (" << g[*e.first].label << ")" << endl;
    //     }
    // }
}

void print_rank_within_cluster(graph& g, int clusterNum){
    struct name_rank {
			double rank;
			string name;
			// 最後のconstを忘れると"instantiated from here"というエラーが出てコンパイルできないので注意
			bool operator<( const name_rank& right ) const {
				return rank < right.rank ;
			}
		};
    vector<name_rank> ranking_list;

    vertex_iterator i,j;
    for (boost::tie(i, j) = vertices(g); *i< xNum ; i++) {
        name_rank a;
        if(g[*i].rx > 0){
            a.rank = g[*i].rx;
            a.name = g[*i].name;
            ranking_list.push_back(a);
        }
    }

    sort(ranking_list.rbegin(), ranking_list.rend());
    cout << "<Cluster = " << clusterNum + 1<< "> (size: " << ranking_list.size() << ")" << endl;
    //top_kがクラスタのサイズより大きかったら変更
    if(top_k > ranking_list.size())top_k = ranking_list.size();
	
    for(int i = 0; i < top_k; i++){
    	cout << i +1 << ": " << ranking_list[i].name << " ... ["<<ranking_list[i].rank  << "]" <<endl;
	}
	cout << endl;
}

void print_cluster(graph& g ){
    vector<vector<string>> cluster(K);
    vertex_iterator i, j;
    for (boost::tie(i, j) = vertices(g); *i< xNum ; i++) {
        cluster[g[*i].belongs_to_cluster].push_back(g[*i].name);
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
