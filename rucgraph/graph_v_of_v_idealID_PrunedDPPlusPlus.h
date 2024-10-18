#pragma once

/*
Li, Rong-Hua, Lu Qin, Jeffrey Xu Yu, and Rui Mao. "Efficient and progressive group steiner tree search." In Proceedings of the 2016 International Conference on Management of Data, pp. 91-106. 2016.
*/

#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp>

#include "graph_hash_of_mixed_weighted_MST_postprocessing.h"
#include <graph_v_of_v_idealID/common_algorithms/graph_v_of_v_idealID_shortest_paths.h>
#include "graph_v_of_v_idealID_save_for_GSTP.h"

#pragma region
pair<vector<int>, vector<int>> graph_v_of_v_idealID_PrunedDPPlusPlus_find_SPs_to_g(graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph, int g_vertex)
{

	/*time complexity: O(|E|+|V|log|V|)*/

	int N = input_graph.size();

	/*add dummy vertex and edges; time complexity: O(|V|)*/
	input_graph.resize(N + 1); // add a dummy vertex N, which is considered as g_vertex
	auto ite = group_graph[g_vertex].end();
	for (auto it = group_graph[g_vertex].begin(); it != ite; it++)
	{

		graph_v_of_v_idealID_add_edge(input_graph, N, it->first, 0); // add dummy edge
	}

	/*time complexity: O(|E|+|V|log|V|)*/
	vector<int> distances; // total vertex and edge weights of paths
	vector<int> predecessors;
	graph_v_of_v_idealID_shortest_paths(input_graph, N, distances, predecessors);

	for (auto it = group_graph[g_vertex].begin(); it != ite; it++)
	{
		graph_hash_of_mixed_weighted_binary_operations_erase(input_graph[it->first], N);
	}
	input_graph.resize(N); // remove dummy vertex N

	distances.resize(N);	// remove dummy vertex N
	predecessors.resize(N); // remove dummy vertex N
	for (int i = 0; i < N; i++)
	{
		if (predecessors[i] == N)
		{
			predecessors[i] = i; // since N is not in predecessors, predecessors[i] points to i, i.e., the path ends at i
		}
	}
	for (size_t i = 0; i < distances.size(); i++)
	{
		if (distances[i] < 0)
		{
			cout << "distance0 " << i << endl;
			exit(1);
		}
	}

	return {predecessors, distances};
}

std::unordered_map<int, pair<vector<int>, vector<int>>> graph_v_of_v_idealID_PrunedDPPlusPlus_find_SPs(graph_v_of_v_idealID &input_graph, graph_v_of_v_idealID &group_graph, std::unordered_set<int> &cumpulsory_group_vertices)
{

	/*time complexity: O(|T||E|+|T||V|log|V|)*/
	std::unordered_map<int, pair<vector<int>, vector<int>>> SPs_to_groups;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		SPs_to_groups[*it] = graph_v_of_v_idealID_PrunedDPPlusPlus_find_SPs_to_g(group_graph, input_graph, *it);
	}

	return SPs_to_groups;
}

#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_find_SPs

#pragma region
struct graph_v_of_v_idealID_PrunedDPPlusPlus_min_node
{
	int v;
	int p; // group_set_ID
	int priority_value;
};
bool operator<(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node const &x, graph_v_of_v_idealID_PrunedDPPlusPlus_min_node const &y)
{
	return x.priority_value > y.priority_value; // < is the max-heap; > is the mean heap; PriorityQueue is expected to be a max-heap of integer values
}
typedef typename boost::heap::fibonacci_heap<graph_v_of_v_idealID_PrunedDPPlusPlus_min_node>::handle_type handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node;
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus priority queue

#pragma region
class graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node
{
	/*this is like the tree T(v,p) in the DPBF paper*/

public:
	int type; // =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge

	int cost; // cost of this tree T(v,p);

	int u; // if this tree is built by grown, then it's built by growing edge (v,u);

	int p1, p2; // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);
};
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node

#pragma region
vector<vector<int>> graph_v_of_v_idealID_PrunedDPPlusPlus_non_overlapped_group_sets(int group_sets_ID_range)
{

	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID;

	time complexity: O(4^|Gamma|), since group_sets_ID_range=2^|Gamma|;

	the original DPBF code use the same method in this function, and thus has the same O(4^|Gamma|) complexity;*/

	vector<vector<int>> non_overlapped_group_sets_IDs(group_sets_ID_range + 1); // <set_ID, non_overlapped_group_sets_IDs>

	for (int i = 1; i <= group_sets_ID_range; i++)
	{ // i is a nonempty group_set ID
		non_overlapped_group_sets_IDs[i] = {};
		for (int j = 1; j < group_sets_ID_range; j++)
		{ // j is another nonempty group_set ID
			if ((i & j) == 0)
			{ // i and j are non-overlapping group sets
				/* The & (bitwise AND) in C or C++ takes two numbers as operands and does AND on every bit of two numbers. The result of AND for each bit is 1 only if both bits are 1.
				https://www.programiz.com/cpp-programming/bitwise-operators */
				non_overlapped_group_sets_IDs[i].push_back(j);
			}
		}
	}

	return non_overlapped_group_sets_IDs;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_non_overlapped_group_sets

#pragma region
vector<vector<int>> graph_v_of_v_idealID_PrunedDPPlusPlus_covered_uncovered_groups(int group_sets_ID_range, std::unordered_set<int> &cumpulsory_group_vertices)
{

	/*time complexity: O(|Gamma|*2^|Gamma|); for each p \in [1,group_sets_ID_range], this function calculate the groups that have not been coverred by p*/

	vector<vector<int>> uncovered_groups(group_sets_ID_range + 1);

	for (int p = 1; p <= group_sets_ID_range; p++)
	{

		vector<int> unc_groups;

		int pow_num = 0;
		for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
		{
			int id = pow(2, pow_num);
			if ((id | p) != p)
			{							   // id is not covered by p
				unc_groups.push_back(*it); // *it is a group not covered by p
			}
			pow_num++;
		}
		uncovered_groups[p] = unc_groups;
	}

	return uncovered_groups;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_covered_uncovered_groups

#pragma region
std::unordered_map<int, std::unordered_map<int, int>> graph_v_of_v_idealID_PrunedDPPlusPlus_virtual_node_distances(graph_v_of_v_idealID &group_graph, std::unordered_set<int> &cumpulsory_group_vertices,
																												   std::unordered_map<int, pair<vector<int>, vector<int>>> &SPs_to_groups)
{

	std::unordered_map<int, std::unordered_map<int, int>> virtual_node_distances;

	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		int g1 = *it;
		auto xx = SPs_to_groups.find(g1);
		for (auto it2 = cumpulsory_group_vertices.begin(); it2 != cumpulsory_group_vertices.end(); it2++)
		{
			int g2 = *it2;
			if (g1 <= g2)
			{
				int distance = std::numeric_limits<int>::max();
				auto pointer_begin = group_graph[g2].begin(), pointer_end = group_graph[g2].end();
				for (auto it2 = pointer_begin; it2 != pointer_end; it2++)
				{
					int dis = xx->second.second[it2->first];
					if (dis < distance)
					{
						distance = dis;
					}
				}
				if (distance < 0)
				{
					cout << "diserr " << g1 << " " << g2 << endl;
					exit(1);
				}

				virtual_node_distances[g1][g2] = distance;
				virtual_node_distances[g2][g1] = distance;
			}
		}
	}

	return virtual_node_distances;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_virtual_node_distances

#pragma region
int graph_v_of_v_idealID_PrunedDPPlusPlus_vertex_group_set_ID(int vertex, graph_v_of_v_idealID &group_graph,
															  std::unordered_set<int> &cumpulsory_group_vertices)
{

	/*time complexity: O(|Gamma|); this function returns the maximum group set ID for a single vertex*/

	int ID = 0;
	int pow_num = 0;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		// cout<<" vid "<<vertex<<" "<<ID<<endl;
		if (graph_v_of_v_idealID_contain_edge(group_graph, *it, vertex))
		{ // vertex is in group *it
			ID = ID + pow(2, pow_num);
		}
		pow_num++;
	}
	return ID;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_vertex_group_set_ID

#pragma region
int graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(int group_vertex,
																								   std::unordered_set<int> &cumpulsory_group_vertices)
{

	/*time complexity: O(|Gamma|); this function returns the maximum group set ID for a single vertex*/

	int pow_num = 0;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		if (*it == group_vertex)
		{
			return pow(2, pow_num);
		}
		pow_num++;
	}
	return 0;
}

#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID

#pragma region
struct graph_v_of_v_idealID_AllPaths_min_node
{
	int v_i, v_j;
	int Xslash;			// group_set_ID
	int priority_value; // W(v_i,v_j,X)
};
bool operator<(graph_v_of_v_idealID_AllPaths_min_node const &x, graph_v_of_v_idealID_AllPaths_min_node const &y)
{
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap; PriorityQueue is expected to be a max-heap of integer values
}
typedef typename boost::heap::fibonacci_heap<graph_v_of_v_idealID_AllPaths_min_node>::handle_type handle_graph_v_of_v_idealID_AllPaths_min_node;

std::unordered_map<string, int> graph_v_of_v_idealID_AllPaths(std::unordered_set<int> &cumpulsory_group_vertices, std::unordered_map<int, std::unordered_map<int, int>> &virtual_node_distances)
{

	std::unordered_map<string, int> W; // String is ("v_i" + "_" + "v_j" + "_" + "group set ID")

	boost::heap::fibonacci_heap<graph_v_of_v_idealID_AllPaths_min_node> Q; // min queue
	std::unordered_map<string, int> Q_priorities;
	std::unordered_map<string, handle_graph_v_of_v_idealID_AllPaths_min_node> Q_handles; // key is String is ("v_i" + "_" + "v_j" + "_" + "group set ID")

	/*D records the popped out optimal subtrees; String is ("v_i" + "_" + "v_j" + "_" + "group set ID") */
	std::unordered_set<string> D;

	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		int p = *it;
		// int Xslash = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(p, cumpulsory_group_vertices);
		graph_v_of_v_idealID_AllPaths_min_node x;

		// x.v_i = p;
		// x.v_j = p;
		// x.Xslash = Xslash;
		// x.priority_value = 0;
		// string handle_ID = to_string(p) + "_" + to_string(p) + "_" + to_string(Xslash);
		// Q_handles[handle_ID] = Q.push(x);
		// Q_priorities[handle_ID] = 0;

		/*the following code for Xslash=0 is not in 2016 paper, but is necessary for computing every W(v_i,v_j,X)*/
		x.v_i = p;
		x.v_j = p;
		x.Xslash = 0;
		x.priority_value = 0;
		string handle_ID = to_string(p) + "_" + to_string(p) + "_" + to_string(0);
		Q_handles[handle_ID] = Q.push(x);
		Q_priorities[handle_ID] = 0;
	}

	while (Q.size() > 0)
	{

		graph_v_of_v_idealID_AllPaths_min_node top_node = Q.top();
		int v_i = top_node.v_i, v_j = top_node.v_j, Xslash = top_node.Xslash;
		int cost = top_node.priority_value;
		Q.pop();

		string handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(Xslash);
		W[handle_ID] = cost;
		D.insert(handle_ID);
		Q_handles.erase(handle_ID);
		Q_priorities.erase(handle_ID);
		// cout << "Q pop " + handle_ID << " priority " << cost << endl;

		/*the following code is not in 2016 paper, but is necessary for computing every W(v_i,v_j,X)*/
		for (int i = 0; i < Xslash; i++)
		{
			if ((i | Xslash) == Xslash)
			{
				handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(i);
				if (W.count(handle_ID) == 0)
				{
					W[handle_ID] = cost;
				}
				else
				{
					if (W[handle_ID] > cost)
					{
						W[handle_ID] = cost;
					}
				}
			}
		}

		for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
		{
			int p = *it;
			int p_setID = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(p, cumpulsory_group_vertices);

			if ((p_setID | Xslash) != Xslash)
			{ // p_setID is not covered by Xslash

				int Xwave = Xslash + p_setID;

				int cost_wave = cost + virtual_node_distances[v_j][p];
				handle_ID = to_string(v_i) + "_" + to_string(p) + "_" + to_string(Xwave);

				if (D.count(handle_ID) > 0)
				{
					continue;
				}

				graph_v_of_v_idealID_AllPaths_min_node x;
				x.v_i = v_i;
				x.v_j = p;
				x.Xslash = Xwave;
				x.priority_value = cost_wave;

				if (Q_handles.count(handle_ID) == 0)
				{
					Q_handles[handle_ID] = Q.push(x);
					Q_priorities[handle_ID] = cost_wave;
					// cout << "Q push " + handle_ID << " priority " << cost_wave << endl;
				}
				else
				{
					if (cost_wave < Q_priorities[handle_ID])
					{
						Q.update(Q_handles[handle_ID], x); // O(1) for decrease key
						Q_priorities[handle_ID] = cost_wave;
						// cout << "Q update " + handle_ID << " priority " << cost_wave << endl;
					}
				}
			}
		}
	}

	// cout << "W.size(): " << W.size() << endl;

	return W;
}
#pragma endregion graph_v_of_v_idealID_AllPaths

#pragma region
graph_hash_of_mixed_weighted graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(int root_v, int root_p, graph_v_of_v_idealID &input_graph,
																			  vector<std::unordered_map<int, graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node>> &trees)
{

	/*this function builds tree T(v,p) at a cost of O(|V|)*/

	graph_hash_of_mixed_weighted solution_tree;

	std::queue<pair<int, int>> waited_to_processed_trees; // <v, p>
	waited_to_processed_trees.push({root_v, root_p});

	while (waited_to_processed_trees.size() > 0)
	{

		int v = waited_to_processed_trees.front().first, p = waited_to_processed_trees.front().second;
		waited_to_processed_trees.pop();

		/*insert v*/
		graph_hash_of_mixed_weighted_add_vertex(solution_tree, v, 0);

		auto pointer_trees_v_p = trees[v].find(p);
		int form_type = pointer_trees_v_p->second.type;
		if (form_type == 0)
		{ // T(v,p) is a single vertex
		}
		else if (form_type == 1)
		{ // T(v,p) is formed by grow
			int u = pointer_trees_v_p->second.u;
			waited_to_processed_trees.push({u, p});
			/*insert (u,v); no need to insert weight of u here, which will be inserted later for T(u,p)*/
			int c_uv = graph_v_of_v_idealID_edge_weight(input_graph, u, v);
			graph_hash_of_mixed_weighted_add_edge(solution_tree, u, v, c_uv);
		}
		else
		{ // T(v,p) is formed by merge
			int p1 = pointer_trees_v_p->second.p1, p2 = pointer_trees_v_p->second.p2;
			waited_to_processed_trees.push({v, p1});
			waited_to_processed_trees.push({v, p2});
		}
	}

	return solution_tree;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree

#pragma region
int graph_v_of_v_idealID_PrunedDPPlusPlus_LB_procedure(int v, int X, int cost, int group_sets_ID_range, vector<vector<int>> &uncovered_groups,
													   std::unordered_map<int, pair<vector<int>, vector<int>>> &SPs_to_groups, std::unordered_map<string, int> &W, std::unordered_map<string, int> &W2)
{

	int lb = 0; // lb should be lower bound cost of a feasible solution contains T(v,X)

	if (group_sets_ID_range != X)
	{
		int X_slash = group_sets_ID_range - X; // X_slash \cup X equals group_sets_ID_range

		int lb1 = INT_MAX, lb2 = -1, lb_one_label = -1;
		vector<int> *pointer_1 = &(uncovered_groups[X]);
		for (auto it = (*pointer_1).begin(); it != (*pointer_1).end(); it++)
		{
			int i = *it;
			int dis_v_i = SPs_to_groups[i].second[v];
			int xxx = INT_MAX;
			for (auto it2 = (*pointer_1).begin(); it2 != (*pointer_1).end(); it2++)
			{
				int j = *it2;
				int dis_v_j = SPs_to_groups[j].second[v];
				if (xxx > dis_v_j)
				{
					xxx = dis_v_j;
				}
				int lb1_value = (dis_v_i + W[to_string(i) + "_" + to_string(j) + "_" + to_string(X_slash)] + dis_v_j) / 2;
				if (lb1 > lb1_value)
				{
					lb1 = lb1_value;
				}
			}
			int lb2_value = (dis_v_i + W2[to_string(i) + "_" + to_string(X_slash)] + xxx) / 2;
			if (lb2 < lb2_value)
			{
				lb2 = lb2_value;
			}
			if (lb_one_label < dis_v_i)
			{
				lb_one_label = dis_v_i;
			}
		}

		if (lb1 < lb2)
		{
			lb = lb2;
		}
		else
		{
			lb = lb1;
		}
		if (lb < lb_one_label)
		{
			lb = lb_one_label;
		}

		// cout << "lb_one_label=" << lb_one_label << endl;
	}

	return cost + lb;
}
#pragma endregion graph_v_of_v_idealID_PrunedDPPlusPlus_LB_procedure

graph_hash_of_mixed_weighted graph_v_of_v_idealID_PrunedDPPlusPlus(graph_v_of_v_idealID &input_graph, graph_v_of_v_idealID &group_graph, std::unordered_set<int> &cumpulsory_group_vertices, int maximum_return_app_ratio, int &RAM_MB, int *pointer2)
{

	/*this function returns the first found feasible solution that has an approximation ratio not larger than maximum_return_app_ratio*/
	int bit_num = 0;

	if (cumpulsory_group_vertices.size() >= 15)
	{
		std::cout << "cumpulsory_group_vertices.size() is too large!" << std::endl;
		exit(1);
	}

	int N = input_graph.size();
	int error_safeguard = 1e-5;
	int inf = std::numeric_limits<int>::max(); // cost of empty tree is inf

	/*finding lowest-weighted paths from groups to vertices; time complexity: O(|T||E|+|T||V|log|V|); return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, pair<vector<int>, vector<int>>> SPs_to_groups = graph_v_of_v_idealID_PrunedDPPlusPlus_find_SPs(input_graph, group_graph, cumpulsory_group_vertices);
	for (auto it = SPs_to_groups.begin(); it != SPs_to_groups.end(); it++)
	{
		bit_num += ((*it).second.first.size() + (*it).second.second.size()) * 4 + sizeof(pair<vector<int>, vector<int>>) + 4 + 3 * 8;
	}

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_v_of_v_idealID_PrunedDPPlusPlus_min_node> Q_T;			   // min queues of trees
	std::unordered_map<string, handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node> Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T
	std::unordered_map<string, int> Q_T_priorities;												   // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_priorities only contain keys that are in Q_T

	/* this is the cost of the best found solution yet */
	graph_hash_of_mixed_weighted best_solu;
	int best_cost = inf;
	//cout << "init " << best_cost << endl;
	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1;																	  // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	vector<vector<int>> non_overlapped_group_sets_IDs = graph_v_of_v_idealID_PrunedDPPlusPlus_non_overlapped_group_sets(group_sets_ID_range); // time complexity: O(4^|Gamma|)
	int xxm = sizeof(vector<int>);
	auto ite = non_overlapped_group_sets_IDs.end();
	for (auto it = non_overlapped_group_sets_IDs.begin(); it != ite; it++)
	{
		bit_num += xxm + (*it).size() * 4;
	}

	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	vector<vector<int>> uncovered_groups = graph_v_of_v_idealID_PrunedDPPlusPlus_covered_uncovered_groups(group_sets_ID_range, cumpulsory_group_vertices);
	auto ite1 = uncovered_groups.end();
	for (auto it = uncovered_groups.begin(); it != ite1; it++)
	{
		bit_num += xxm + (*it).size() * 4;
	}

	std::unordered_map<int, std::unordered_map<int, int>> virtual_distances = graph_v_of_v_idealID_PrunedDPPlusPlus_virtual_node_distances(group_graph, cumpulsory_group_vertices, SPs_to_groups);
	for (auto i = SPs_to_groups.begin(); i != SPs_to_groups.end(); i++)
	{
		/* code */
	}

	xxm = sizeof(std::unordered_map<int, int>);
	auto ite2 = virtual_distances.end();
	for (auto it = virtual_distances.begin(); it != ite2; it++)
	{
		bit_num += xxm + (it->second).size() * 8;
	}

	std::unordered_map<string, int> W = graph_v_of_v_idealID_AllPaths(cumpulsory_group_vertices, virtual_distances); // String is ("v_i" + "_" + "v_j" + "_" + "group set ID")
	std::unordered_map<string, int> W2;																				 // String is ("v_i" + "_" + "group set ID")
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
	{
		int v_i = *it;
		for (int Xslash = 1; Xslash <= group_sets_ID_range; Xslash++)
		{
			int dis = INT_MAX;
			for (auto it2 = cumpulsory_group_vertices.begin(); it2 != cumpulsory_group_vertices.end(); it2++)
			{
				int v_j = *it2;
				string handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(Xslash);
				if (dis > W[handle_ID])
				{
					dis = W[handle_ID];
				}
			}
			string handle_ID = to_string(v_i) + "_" + to_string(Xslash);
			W2[handle_ID] = dis;
		}
	}
	bit_num += 2 * 8 + W.size() * 12;
	bit_num += 2 * 8 + W2.size() * 12;

	/*initialize trees with vertices; time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees*/
	vector<std::unordered_map<int, graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node>> trees(N); // <v, <p, T(v,p)>>; this cannot be changed to vector of vectors
	xxm = sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node);
	for (int v = 0; v < N; v++)
	{
		int group_set_ID_v = graph_v_of_v_idealID_PrunedDPPlusPlus_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/

		for (int p = 1; p <= group_set_ID_v; p++)
		{ // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v)
			{ // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
				Q_T_priorities[handle_ID] = x.priority_value;
				// cout << "initial Q push " + handle_ID << " priority: " << x.priority_value << endl;

				bit_num += xxm + 2 * 8;
			}
		}
	}

	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;

	// cout << "group_sets_ID_range:" << group_sets_ID_range << endl;
	int show = 1;
	int Q_T_max_size = 0;
	graph_v_of_v_idealID_PrunedDPPlusPlus_min_node top_node;
	/*Big while loop*/
	int counts = 0;
	while (Q_T.size() > 0)
	{ // at most 2^|Gamma|*V loops
		counts++;

		int Q_T_size = Q_T.size();
		if (Q_T_size > Q_T_max_size)
		{
			Q_T_max_size = Q_T_size;
		}

		top_node = Q_T.top();
		int v = top_node.v, X = top_node.p;
		

		int v_X_tree_cost = trees[v][X].cost;

		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements
				   //	cout<<"V:"<<v<<" X: "<<X<<endl;
		string handle_ID = to_string(v) + "_" + to_string(X);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T
		Q_T_priorities.erase(handle_ID);

		if (X == group_sets_ID_range)
		{ // T(v,p) covers all groups
			// cout << "X:" << X << endl;
			bit_num += Q_T_max_size * (sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8 + sizeof(handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8); // Q_T + Q_T_handles
			RAM_MB = bit_num / 1024 / 1024;

			graph_hash_of_mixed_weighted feasible_solu = graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(v, X, input_graph, trees);
			return feasible_solu;
		}

		// cout << "Q pop " + handle_ID << endl;
		D.insert(handle_ID); // optimal T(v,p) has been found

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		int feasible_solu_cost = v_X_tree_cost;
		auto ite3 = uncovered_groups[X].end();
		
		for (auto it = uncovered_groups[X].begin(); it != ite3; it++)
		{
			if (SPs_to_groups[*it].second[v] == inf)
			{
				feasible_solu_cost = inf;
				break;
			}
			else
			{

				feasible_solu_cost = feasible_solu_cost + SPs_to_groups[*it].second[v];
			}
		}
		
		if (feasible_solu_cost < best_cost)
		{
			graph_hash_of_mixed_weighted feasible_solu = graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(v, X, input_graph, trees);
			for (auto it = uncovered_groups[X].begin(); it != uncovered_groups[X].end(); it++)
			{
				int g_id = *it;

				/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
				int v_start = v;
				int pre = SPs_to_groups[g_id].first[v_start];
				
				// cout<<"gid "<<g_id<<" vstart "<<v_start<<" pre "<<pre<<endl;
				while (v_start != pre)
				{
					int ec = graph_v_of_v_idealID_edge_weight(input_graph, v_start, pre);
					graph_hash_of_mixed_weighted_add_vertex(feasible_solu, pre, 0);
					// cout<<"add cost "<<ec<<endl;
					graph_hash_of_mixed_weighted_add_edge(feasible_solu, v_start, pre, ec);
					v_start = pre;
					pre = SPs_to_groups[g_id].first[v_start];
				}
			}

			best_solu = graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(feasible_solu);
			best_cost = graph_hash_of_mixed_weighted_sum_of_ec(best_solu);
		}


		// cout<<"after one label best "<<best_cost<<endl;
		/*since an optimal solution has not been popped out, the optimal cost must be larger than or equal to v_X_tree_cost*/
		if (v_X_tree_cost > 0)
		{
			// cout << "X:" << X << endl;
			int ratio = best_cost / v_X_tree_cost;
			// cout << "ratio:" << ratio << " v_X_tree_cost:" << v_X_tree_cost << endl;
			if (ratio <= maximum_return_app_ratio + 1e-5)
			{ // this feasible solution can be returned
				// std::cout<<"count "<<counts<<std::endl;
				*pointer2 += counts;
				bit_num += Q_T_max_size * (sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8 + sizeof(handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8); // Q_T + Q_T_handles
				RAM_MB = bit_num / 1024 / 1024;

				return best_solu;
			}
		}

		/*Lines 17-19 in PrunedDP++*/
		int X_slash = group_sets_ID_range - X; // p \cup X_slash equals group_sets_ID_range
		handle_ID = to_string(v) + "_" + to_string(X_slash);
		if (D.count(handle_ID) > 0)
		{

			int merged_tree_cost = v_X_tree_cost + trees[v][X_slash].cost;

			int cost_Tvp1_cup_p2 = inf;
			if (trees[v].count(group_sets_ID_range) > 0)
			{
				cost_Tvp1_cup_p2 = trees[v][group_sets_ID_range].cost;
			}
			if (merged_tree_cost < cost_Tvp1_cup_p2)
			{ // update tree T(v,group_sets_ID_range)
				// cout<<"merge"<<endl;
				graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node node_x;
				node_x.cost = merged_tree_cost;
				node_x.type = 2;
				node_x.p1 = X;
				node_x.p2 = X_slash;
				trees[v][group_sets_ID_range] = node_x;

				int lb = merged_tree_cost;

				if (merged_tree_cost > best_cost + error_safeguard)
				{ // cannot have >= here, since best_solu may not be in Q; if >=, then best_solu may never be in Q
					continue;
				}
				if (merged_tree_cost < best_cost)
				{
					best_cost = merged_tree_cost;
					best_solu = graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(v, group_sets_ID_range, input_graph, trees);
				}

				handle_ID = to_string(v) + "_" + to_string(group_sets_ID_range);
				if (Q_T_priorities.count(handle_ID) == 0)
				{
					graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
					x.v = v;
					x.p = group_sets_ID_range;
					x.priority_value = merged_tree_cost;
					Q_T_handles[handle_ID] = Q_T.push(x);
					Q_T_priorities[handle_ID] = x.priority_value;
					// cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
				}
				else
				{
					if (lb < Q_T_priorities[handle_ID])
					{
						graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
						x.v = v;
						x.p = group_sets_ID_range;
						x.priority_value = merged_tree_cost;
						Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
						Q_T_priorities[handle_ID] = x.priority_value;
					}
				}
			}
			continue;
		}

		/* Optimal-Tree Decomposition Theorem (Theorem 1 in 2016 SIGMOD paper)*/
		if (v_X_tree_cost < best_cost / 2 + error_safeguard)
		{

			// cout << "here2" << endl;

			/*grow*/
			auto ite4 = input_graph[v].end();
			for (auto it = input_graph[v].begin(); it != ite4; it++)
			{

				/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
				int u = it->first;
				int cost_euv = it->second;
				int grow_tree_cost = trees[v][X].cost + cost_euv;

				handle_ID = to_string(u) + "_" + to_string(X);
				// cout << "grow " << handle_ID <<  " grow_tree_cost " << grow_tree_cost << " best_cost " << best_cost << endl;
				if (D.count(handle_ID) > 0)
				{
					continue;
				}

				int T_up_cost;
				if (trees[u].count(X) == 0)
				{
					T_up_cost = inf;
				}
				else
				{
					T_up_cost = trees[u][X].cost;
				}
				if (grow_tree_cost < T_up_cost)
				{
					graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node node_x;
					node_x.cost = grow_tree_cost;
					node_x.type = 1;
					node_x.u = v;
					trees[u][X] = node_x;

					int lb = graph_v_of_v_idealID_PrunedDPPlusPlus_LB_procedure(u, X, grow_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);

					// cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;

					if (lb > best_cost + error_safeguard)
					{
						// cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;
						continue;
					}

					if (Q_T_priorities.count(handle_ID) == 0)
					{
						graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
						x.v = u;
						x.p = X;
						x.priority_value = grow_tree_cost;
						Q_T_handles[handle_ID] = Q_T.push(x);
						Q_T_priorities[handle_ID] = x.priority_value;

						// cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
					}
					else
					{
						if (grow_tree_cost < Q_T_priorities[handle_ID])
						{
							graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
							x.v = u;
							x.p = X;
							x.priority_value = grow_tree_cost;
							Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
							Q_T_priorities[handle_ID] = x.priority_value;
						}
					}
				}
			}

			// cout << "here3" << endl;

			/*merge*/
			int p1 = X;
			auto ite5 = non_overlapped_group_sets_IDs[p1].end();
			for (auto it = non_overlapped_group_sets_IDs[p1].begin(); it != ite5; it++)
			{
				// cout<<"merge2"<<endl;
				int p2 = *it; // p2 is not overlapped with p1
				if (D.count(to_string(v) + "_" + to_string(p2)) > 0)
				{ // only merge optimal T(v,p2), this is recorded in the pseudo code of 2016 paper

					int p1_cup_p2 = p1 + p2;
					handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);

					if (D.count(handle_ID) > 0)
					{
						continue;
					}

					int merged_tree_cost = v_X_tree_cost + trees[v][p2].cost;

					// cout << "merge " << handle_ID << " merged_tree_cost " << merged_tree_cost << " best_cost " << best_cost << endl;

					int cost_Tvp1_cup_p2;
					if (trees[v].count(p1_cup_p2) == 0)
					{
						cost_Tvp1_cup_p2 = inf;
					}
					else
					{
						cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
					}

					if (merged_tree_cost < cost_Tvp1_cup_p2)
					{ // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)
						graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node node_x;
						node_x.cost = merged_tree_cost;
						node_x.type = 2;
						node_x.p1 = p1;
						node_x.p2 = p2;
						trees[v][p1_cup_p2] = node_x;

						/*Conditional Tree Merging Theorem (Theorem 2 in 2016 SIGMOD paper);
						note the following place of error_safeguard; this if condition should be met even if error_safeguard=inf*/
						if (merged_tree_cost <= best_cost * 2 / 3 + error_safeguard)
						{ // error_safeguard is error

							int lb = graph_v_of_v_idealID_PrunedDPPlusPlus_LB_procedure(v, p1_cup_p2, merged_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);
							// cout << "merge cost_Tvp1_cup_p2=" << cost_Tvp1_cup_p2 << " lb_slash=" << lb_slash << " lb=" << lb << " best_cost=" << best_cost << endl;

							if (lb > best_cost + error_safeguard)
							{
								continue;
							}
							if (p1_cup_p2 == group_sets_ID_range && merged_tree_cost < best_cost)
							{
								best_cost = merged_tree_cost;
								best_solu = graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(v, group_sets_ID_range, input_graph, trees);
							}

							if (Q_T_priorities.count(handle_ID) == 0)
							{
								graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
								x.v = v;
								x.p = p1_cup_p2;
								x.priority_value = merged_tree_cost;
								Q_T_handles[handle_ID] = Q_T.push(x);
								Q_T_priorities[handle_ID] = x.priority_value;
								// cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
							}
							else
							{
								if (merged_tree_cost < Q_T_priorities[handle_ID])
								{ // if priority is lb values, then here is lv < < Q_T_priorities[handle_ID]
									graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
									x.v = v;
									x.p = p1_cup_p2;
									x.priority_value = merged_tree_cost;
									Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
									Q_T_priorities[handle_ID] = x.priority_value;
									// cout << "Q update " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
								}
							}
						}
					}
				}
			}
		}
	}

	std::cout << "inside graph_v_of_v_idealID_PrunedDPPlusPlus did not find a feasible solution!" << std::endl;
	graph_v_of_v_idealID_save_for_GSTP("PrunedDPPlusPlus_irregular.txt", input_graph, group_graph, cumpulsory_group_vertices);
	getchar();
	exit(1);
}