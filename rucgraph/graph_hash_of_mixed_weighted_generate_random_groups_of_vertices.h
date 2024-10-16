#pragma once

#include <boost/random.hpp>
using namespace std;
void graph_hash_of_mixed_weighted_generate_random_groups_of_vertices(int G, int g_size_min, int g_size_max,
																	 graph_hash_of_mixed_weighted &input_graph, int group_vertex_start_ID,
																	 std::unordered_set<int> &generated_group_vertices, graph_hash_of_mixed_weighted &generated_group_graph)
{

	/*time complexity: O(|G||V|)*/

	std::time_t now = std::time(0);
	boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
	boost::random::uniform_int_distribution<> dist{g_size_min, g_size_max};

	/*time complexity: O(|V|)*/
	std::vector<int> all_v;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++)
	{
		int v = it->first;
		graph_hash_of_mixed_weighted_add_vertex(generated_group_graph, v, 0); // have same v
		all_v.push_back(v);
	}

	/*add groups*/
	/*time complexity: O(|G||V|)*/
	int add_group_num = 0;
	while (add_group_num < G)
	{
		int group_size = dist(gen); // generate int random number
		std::vector<int> to_be_linked_v = all_v;
		std::vector<int> linked_vertices; // enlarge linked size until it become enough;
		while (linked_vertices.size() < group_size)
		{
			boost::random::uniform_int_distribution<> dist2{0, (int)(to_be_linked_v.size() - 1)};
			int randID = dist2(gen);											   // choose from to be linked
			linked_vertices.insert(linked_vertices.end(), to_be_linked_v[randID]); // add vertex to group
			to_be_linked_v.erase(to_be_linked_v.begin() + randID);				   // delete vertex from to_be_linked
		}

		// add this group
		int group_vertex = group_vertex_start_ID + add_group_num; // group insert after vertexs;
		generated_group_vertices.insert(group_vertex);
		for (int j = 0; j < linked_vertices.size(); j++)
		{
			graph_hash_of_mixed_weighted_add_edge(generated_group_graph, linked_vertices[j], group_vertex, 1);
			// edge between group vertex(stands for group) and linked_vertex
		}

		add_group_num++;
	}
}

void graph_hash_of_mixed_weighted_generate_community_groups_of_vertices(int G, int g_size_min, int g_size_max,
																		graph_hash_of_mixed_weighted &input_graph, int group_vertex_start_ID,
																		std::unordered_set<int> &generated_group_vertices, graph_hash_of_mixed_weighted &generated_group_graph, int percent, int *belong, std::vector<std::vector<int>> &community)
{

	/*time complexity: O(|G||V|)*/

	std::time_t now = std::time(0);
	boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
	boost::random::uniform_int_distribution<> dist{g_size_min, g_size_max};

	/*time complexity: O(|V|)*/
	std::vector<int> all_v;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++)
	{
		int v = it->first;
		graph_hash_of_mixed_weighted_add_vertex(generated_group_graph, v, 0); // have same v
		all_v.push_back(v);
	}
	int nums = all_v.size();
	/*add groups*/
	/*time complexity: O(|G||V|)*/
	int add_group_num = 0;
	while (add_group_num < G)
	{

		int group_size = dist(gen); // generate int random number
		std::vector<int> to_be_linked_v = all_v;
		std::vector<int> linked_vertices; // enlarge linked size until it become enough;
		boost::random::uniform_int_distribution<> distf{0, nums - 1};
		int first = distf(gen);
		std::vector<int>c  =community[belong[first]];
		
		while (linked_vertices.size() < group_size * percent / 100&&c.size())
		{
			boost::random::uniform_int_distribution<> dist2{0, (int)(c.size() - 1)};
			int randID = dist2(gen);											   // choose from to be linked
			int need = c[randID];
			linked_vertices.insert(linked_vertices.end(), c[randID]); // add vertex to group
			to_be_linked_v.erase(to_be_linked_v.begin() + need);				   // delete vertex from to_be_linked
			c.erase(c.begin()+randID);
			
		}
		while (linked_vertices.size() < group_size)
		{

			boost::random::uniform_int_distribution<> dist2{0, (int)(to_be_linked_v.size() - 1)};
			int randID = dist2(gen);											   // choose from to be linked
			linked_vertices.insert(linked_vertices.end(), to_be_linked_v[randID]); // add vertex to group
			to_be_linked_v.erase(to_be_linked_v.begin() + randID);				   // delete vertex from to_be_linked
	
		}

		// add this group
		int group_vertex = group_vertex_start_ID + add_group_num; // group insert after vertexs;
		generated_group_vertices.insert(group_vertex);
		for (int j = 0; j < linked_vertices.size(); j++)
		{
			graph_hash_of_mixed_weighted_add_edge(generated_group_graph, linked_vertices[j], group_vertex, 1);
			// edge between group vertex(stands for group) and linked_vertex
		}

		add_group_num++;
	}
}


