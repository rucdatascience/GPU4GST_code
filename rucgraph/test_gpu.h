#pragma once



#include <chrono>
#include <queue>
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


bool this_is_a_feasible_solution_gpu(graph_hash_of_mixed_weighted& solu, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& group_vertices) {

	/*time complexity O(|V_solu|+|E_solu|)*/
	if (graph_hash_of_mixed_weighted_connected_components(solu).size() != 1) { // it's not connected
		cout << "this_is_a_feasible_solution: solu is disconnected!" << endl;
		return false;
	}

	for (auto it = group_vertices.begin(); it != group_vertices.end(); it++) {
		int g = *it;
		bool covered = false;
		for (auto it2 = solu.hash_of_vectors.begin(); it2 != solu.hash_of_vectors.end(); it2++) {
			int v = it2->first;
			if (graph_hash_of_mixed_weighted_contain_edge(group_graph, v, g)) {
				covered = true;
				break;
			}
		}
		if (covered == false) {
			cout << "this_is_a_feasible_solution: a group is not covered!" << endl;
			return false;
		}
	}

	return true;

}


void test_graph_v_of_v_idealID_DPBF_only_ec_gpu() {

	/*parameters*/
	int iteration_times = 1;
	int V = 10000, E = 50000, G = 4, g_size_min = 1, g_size_max = 4, precision = 0;
	int ec_min = 1, ec_max = 4; // PrunedDP does not allow zero edge weight



	int solution_cost_DPBF_sum = 0, solution_cost_PrunedDPPlusPlus_sum = 0;

	double time_DPBF_avg = 0, time_PrunedDPPlusPlus_avg = 0;

	/*iteration*/
	for (int i = 0; i < iteration_times; i++) {

		cout << "iteration " << i << endl;

		/*input and output*/
		int generate_new_graph = 1;
		int lambda = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);

			graph_hash_of_mixed_weighted_generate_random_groups_of_vertices(G, g_size_min, g_size_max,
				instance_graph, instance_graph.hash_of_vectors.size(), generated_group_vertices, generated_group_graph);//
			
			
			graph_hash_of_mixed_weighted_save_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}
		else {
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}
		cout<<"generate complete"<<endl;
		unordered_map<int, int> vertexID_old_to_new;
		for (int mm = 0; mm < V; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_instance_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(instance_graph, vertexID_old_to_new);
		for (int mm = V; mm < V + G; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_generated_group_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(generated_group_graph, vertexID_old_to_new);
		CSR_graph csr_graph = toCSR(v_instance_graph);
		cout<<"E:"<<csr_graph.E_all<<" v:"<<csr_graph.V<<endl;
		/*graph_v_of_v_idealID_DPBF_only_ec*/
		if (1) {
			int RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = DPBF_GPU(csr_graph,generated_group_vertices,v_generated_group_graph,v_instance_graph);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_DPBF_avg += (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			int cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			cout<<"cost: "<<cost<<endl;
			solution_cost_DPBF_sum += cost;

			if (!this_is_a_feasible_solution_gpu(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}


		/*graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted*/
		if (1) {
			int RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_v_of_v_idealID_PrunedDPPlusPlus(v_instance_graph, v_generated_group_graph, generated_group_vertices, 1, RAM);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_PrunedDPPlusPlus_avg = time_PrunedDPPlusPlus_avg + (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			int cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_PrunedDPPlusPlus_sum = solution_cost_PrunedDPPlusPlus_sum + cost;

			if (!this_is_a_feasible_solution_gpu(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}

		if (solution_cost_DPBF_sum + 1e-8 < solution_cost_PrunedDPPlusPlus_sum) {
			cout << "solution_cost_DPQ_GPU_sum=" << solution_cost_DPBF_sum << endl;
			cout << "solution_cost_PrunedDPPlusPlus_sum=" << solution_cost_PrunedDPPlusPlus_sum << endl;
			cout << "wrong answer" << endl;
			getchar();
		}


	}

	cout << "solution_cost_DPBF_sum=" << solution_cost_DPBF_sum << endl;
	cout << "solution_cost_PrunedDPPlusPlus_sum=" << solution_cost_PrunedDPPlusPlus_sum << endl;

	cout << "time_DPBF_avg=" << time_DPBF_avg << "s" << endl;
	cout << "time_PrunedDPPlusPlus_avg=" << time_PrunedDPPlusPlus_avg << "s" << endl;
}