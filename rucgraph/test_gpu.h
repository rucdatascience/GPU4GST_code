#pragma once

#include <chrono>
#include <queue>
#include <omp.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_connected_components.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted_generate_random_groups_of_vertices.h>
#include <graph_hash_of_mixed_weighted_save_for_GSTP.h>
#include <graph_hash_of_mixed_weighted_read_for_GSTP.h>
#include "graph_hash_of_mixed_weighted_sum_of_nw_ec.h"
#include "graph_v_of_v_idealID_DPBF_only_ec.h"
#include "graph_v_of_v_idealID_PrunedDPPlusPlus.h"
#include "DPQ.cuh"
#include <future>
non_overlapped_group_sets graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(int group_sets_ID_range)
{
	non_overlapped_group_sets s;
	s.length = 0;
	s.non_overlapped_group_sets_IDs_pointer_host.resize(group_sets_ID_range + 3);
	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID;
	time complexity: O(4^|Gamma|), since group_sets_ID_range=2^|Gamma|;
	the original DPBF code use the same method in this function, and thus has the same O(4^|Gamma|) complexity;*/
	// <set_ID, non_overlapped_group_sets_IDs>
	for (int i = 1; i <= group_sets_ID_range; i++)
	{ // i is a nonempty group_set ID
		s.non_overlapped_group_sets_IDs_pointer_host[i] = s.length;
		for (int j = 1; j < group_sets_ID_range; j++)
		{ // j is another nonempty group_set ID
			if ((i & j) == 0)
			{ // i and j are non-overlapping group sets
				/* The & (bitwise AND) in C or C++ takes two numbers as operands and does AND on every bit of two numbers. The result of AND for each bit is 1 only if both bits are 1.
				https://www.programiz.com/cpp-programming/bitwise-operators */
				s.non_overlapped_group_sets_IDs.push_back(j);
				s.length++;
			}
		}
	}
	s.non_overlapped_group_sets_IDs_pointer_host[group_sets_ID_range + 1] = s.length;
	return s;
}
bool this_is_a_feasible_solution_gpu(graph_hash_of_mixed_weighted &solu,  graph_v_of_v_idealID &group_graph,
									 std::unordered_set<int> &group_vertices)
{

	/*time complexity O(|V_solu|+|E_solu|)*/
	if (graph_hash_of_mixed_weighted_connected_components(solu).size() != 1)
	{ // it's not connected
		cout << "this_is_a_feasible_solution: solu is disconnected!" << endl;
		return false;
	}

	for (auto it = group_vertices.begin(); it != group_vertices.end(); it++)
	{
		int g = *it;
		bool covered = false;
		for (auto it2 = solu.hash_of_vectors.begin(); it2 != solu.hash_of_vectors.end(); it2++)
		{
			int v = it2->first;
			if ( graph_v_of_v_idealID_contain_edge(group_graph, g, v))
			{
				covered = true;
				break;
			}
		}
		if (covered == false)
		{
			cout << "this_is_a_feasible_solution: a group is not covered!" << endl;
			return false;
		}
	}

	return true;
}

void test_graph_v_of_v_idealID_DPBF_only_ec_gpu()
{
	cout << "start.. " << endl;
	/*parameters*/
	int iteration_times = 1;
	int V = 300000, E = 800000, G = 6, g_size_min = 1, g_size_max = 4, precision = 0;
	int ec_min = 1, ec_max = 4; // PrunedDP does not allow zero edge weight

	int solution_cost_GPU_sum = 0, solution_cost_PrunedDPPlusPlus_sum = 0, solution_cost_multi_GPU_sum = 0;
	std::vector<std::vector<int>> community;
	double time_GPU_one_avg = 0, time_PrunedDPPlusPlus_avg = 0, time_GPU_multi_avg = 0;
	int p_gpu = 0, p_cpu = 0;
	int *pointer1 = &p_gpu, *pointer2 = &p_cpu;
	int generate_new_graph = 0;
	int lambda = 1;
	std::unordered_set<int> generated_group_vertices;
	graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
	std::vector<std::vector<int>> inquire;
	graph_v_of_v_idealID v_generated_group_graph, v_instance_graph;
	/*if (generate_new_graph == 1)
	{
		instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);

		graph_hash_of_mixed_weighted_generate_random_groups_of_vertices(G, g_size_min, g_size_max,
																		instance_graph, instance_graph.hash_of_vectors.size(), generated_group_vertices, generated_group_graph); //

		graph_hash_of_mixed_weighted_save_for_GSTP("simple_iterative_tests.text", instance_graph,
												   generated_group_graph, generated_group_vertices, lambda);
		graph_hash_of_mixed_weighted_read_for_GSTP_1(data_name+"_in", instance_graph,
												 generated_group_graph, generated_group_vertices, lambda, community);
	cout << "community size " << community.size() << endl;
	graph_hash_of_mixed_weighted_read_for_Group(data_name+".g", instance_graph,
												generated_group_graph, generated_group_vertices);
	cout << "Group graph size " << generated_group_graph.hash_of_vectors.size() << endl;
	graph_hash_of_mixed_weighted_read_for_inquire(data_name+".g5", inquire);
	}*/
	string data_name = "Amazon";
	read_input_graph(data_name + ".in", v_instance_graph, community);
	cout << "read input complete" << endl;
	read_Group(data_name + ".g", v_instance_graph, v_generated_group_graph);

	cout << "read group complete " << v_generated_group_graph.size() << endl;
	read_inquire(data_name + ".g7", inquire);

	iteration_times = inquire.size();
	cout << "inquires size " << inquire.size()<<" G = "<<inquire[0].size() << endl;
	V = v_instance_graph.size();
	CSR_graph csr_graph = toCSR(v_instance_graph);
	cout << "E: " << csr_graph.E_all << " V: " << csr_graph.V << endl;
	int *belong = new int[csr_graph.V];
	for (size_t i = 0; i < community.size(); i++)
	{
		for (size_t j = 0; j < community[i].size(); j++)
		{
			belong[community[i][j]] = i;
		}
	}

	G = inquire[0].size();
	cout<<"gsize "<<G<<endl;
	int group_sets_ID_range = pow(2, G) - 1;
	non_overlapped_group_sets s = graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(group_sets_ID_range);
	/*iteration*/
	cout << "------------------------------------------------------------" << endl;
	int rounds = 0;
	for (int i =88; i < 89; i++)
	{
			rounds++;
		cout << "iteration " << i << endl;

		// generated_group_graph.clear();
		generated_group_vertices.clear();
		for (size_t j = 0; j < inquire[i].size(); j++)
		{	
			generated_group_vertices.insert(inquire[i][j] - V);
		}

		/*graph_hash_of_mixed_weighted_generate_community_groups_of_vertices(G, g_size_min, g_size_max,
																		   instance_graph, instance_graph.hash_of_vectors.size(), generated_group_vertices, generated_group_graph, 20, belong, community); //
*/

		std::cout << "get inquire complete" << std::endl;
		if (0)
		{
			int RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_v_of_v_idealID_PrunedDPPlusPlus(v_instance_graph, v_generated_group_graph, generated_group_vertices, 1, RAM, pointer2);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_PrunedDPPlusPlus_avg += runningtime ;

			// graph_hash_of_mixed_weighted_print(solu);

			int cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_PrunedDPPlusPlus_sum = solution_cost_PrunedDPPlusPlus_sum + cost;
			cout << "PLus PLus " << cost << " time " << runningtime << endl;
					graph_hash_of_mixed_weighted_print(solu);
			if (!this_is_a_feasible_solution_gpu(solu, v_generated_group_graph, generated_group_vertices))
			{
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				//  exit(1);
			}
			cout << "------------------------------------------------------------" << endl;
		}
		if (1)
		{
			// CSR_graph part_graph[3];
			// 	toCSR_three(part_graph, v_instance_graph, community, c_size);
			node **host_tree;
			int c_size[3];
			int height = csr_graph.V, width = group_sets_ID_range + 1;
			host_tree = new node *[height];
			node *host_tree_one_d = new node[height * width];
			for (size_t i = 0; i < height; i++)
			{
				host_tree[i] = &host_tree_one_d[i * width];
			}
			int cost;
			int RAM, *real_cost = &cost;
			auto begin = std::chrono::high_resolution_clock::now();
/*#pragma omp parallel for
			for (int i = 0; i < 3; ++i)
			{
				// DPBF_GPU_part(host_tree, part_graph[i], generated_group_vertices, v_generated_group_graph, v_instance_graph, pointer1, real_cost, community, c_size, i + 1, s);
			}*/
			// double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			// double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
			double runningtime;
			graph_hash_of_mixed_weighted solu = DPBF_GPU(host_tree, host_tree_one_d, csr_graph, generated_group_vertices, v_generated_group_graph, v_instance_graph, pointer1, real_cost, belong, c_size, s, &runningtime);
			auto end = std::chrono::high_resolution_clock::now();
			runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
			time_GPU_one_avg += runningtime;

			// graph_hash_of_mixed_weighted_print(solu);

			cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_GPU_sum += *real_cost;
			cout << "one form tree cost: " << cost << " full time " << runningtime << endl;
//graph_hash_of_mixed_weighted_print(solu);
			if (!this_is_a_feasible_solution_gpu(solu, v_generated_group_graph, generated_group_vertices))
			{
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
			}
			cout << "------------------------------------------------------------" << endl;
		}
		if (1)
		{
			node **host_tree;
			int height = csr_graph.V, width = group_sets_ID_range + 1;
			host_tree = new node *[height];
			node *host_tree_one_d = new node[height * width];
			for (size_t i = 0; i < height; i++)
			{
				host_tree[i] = &host_tree_one_d[i * width];
			}
			int rcost;
			int RAM, *real_cost = &rcost;

			auto begin = std::chrono::high_resolution_clock::now();
			// double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			// double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
			double runningtime;
			// 		  std::future<graph_hash_of_mixed_weighted> result = std::async(DPBF_GPU_T,host_tree, host_tree_one_d, std::ref(csr_graph),std::ref(generated_group_vertices) ,std::ref(v_generated_group_graph) ,std::ref(v_instance_graph) , pointer1, real_cost, belong, std::ref(community), s, &runningtime);
			// graph_hash_of_mixed_weighted product = result.get();

			std::future<graph_hash_of_mixed_weighted> result[4];
			graph_hash_of_mixed_weighted solutions[4];
			// result[0] = std::async(DPBF_GPU_T,host_tree, host_tree_one_d, std::ref(csr_graph),std::ref(generated_group_vertices) ,std::ref(v_generated_group_graph) ,std::ref(v_instance_graph) ,pointer1, real_cost,  belong, std::ref(community), s, &runningtime);

			// for (size_t i = 1; i < 4; i++)
			// {
			// 	result[i] = std::async(DPBF_GPU_Part,host_tree, host_tree_one_d, std::ref(csr_graph),std::ref(generated_group_vertices) ,std::ref(v_generated_group_graph) ,std::ref(v_instance_graph) ,  belong, std::ref(community), s,i);

			// }
			// for (size_t i = 1; i < 4; i++)
			// {
			// 	solutions[i] =  result[i].get();
			// }

			// solutions[0] = result[0].get();
			solutions[0] = DPBF_GPU_T(host_tree, host_tree_one_d, csr_graph, generated_group_vertices, v_generated_group_graph, v_instance_graph, pointer1, real_cost, belong, community, s, &runningtime);

			auto end = std::chrono::high_resolution_clock::now();
		 runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
			time_GPU_multi_avg += runningtime;
			graph_hash_of_mixed_weighted solu;

			int cost = graph_hash_of_mixed_weighted_sum_of_ec(solutions[0]);
			solution_cost_multi_GPU_sum += *real_cost;
			cout << "community form tree cost: " << cost << " full time " << runningtime << endl;

			if (!this_is_a_feasible_solution_gpu(solutions[0], v_generated_group_graph, generated_group_vertices))
			{
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solutions[0]);
			}
			cout << "------------------------------------------------------------" << endl;
		}

	}
	cudaFree(csr_graph.all_edge);
	cudaFree(csr_graph.all_pointer);
	cudaFree(csr_graph.all_edge_weight);
	cout << "solution_cost_CPU  _sum=" << solution_cost_PrunedDPPlusPlus_sum << endl;
	cout << "solution_cost_GPU_1_sum=" << solution_cost_GPU_sum << endl;
	cout << "solution_cost_GPU_2_sum=" << solution_cost_multi_GPU_sum << endl;
	cout << "time_CPU  _avg=" << time_PrunedDPPlusPlus_avg / rounds<< "s" << endl;
	cout << "time_GPU_1_avg=" << time_GPU_one_avg / rounds<< "s" << endl;
	cout << "time_GPU_2_avg=" << time_GPU_multi_avg/ rounds << "s" << endl;

}