#pragma once

/*
this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.

time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )
*/

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <build_in_progress/GST/test.h>


int main()
{
	test_graph_v_of_v_idealID_DPBF_only_ec();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/

#include<test.cuh>


bool this_is_a_feasible_solution(graph_hash_of_mixed_weighted& solu, graph_hash_of_mixed_weighted& group_graph,
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


void test_graph_v_of_v_idealID_DPBF_only_ec() {

	/*parameters*/
	int iteration_times = 100;
	int V = 1000, E = 5000, G = 5, g_size_min = 10, g_size_max = 50, precision = 3;
	double ec_min = 0.001, ec_max = 1; // PrunedDP does not allow zero edge weight



	double solution_cost_DPBF_sum = 0, solution_cost_PrunedDPPlusPlus_sum = 0;

	double time_DPBF_avg = 0, time_PrunedDPPlusPlus_avg = 0;

	/*iteration*/
	for (int i = 0; i < iteration_times; i++) {

		cout << "iteration " << i << endl;

		/*input and output*/
		int generate_new_graph = 1;
		double lambda = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_generate_random_groups_of_vertices(G, g_size_min, g_size_max,
				instance_graph, instance_graph.hash_of_vectors.size(), generated_group_vertices, generated_group_graph);
			graph_hash_of_mixed_weighted_save_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}
		else {
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}

		unordered_map<int, int> vertexID_old_to_new;
		for (int mm = 0; mm < V; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_instance_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(instance_graph, vertexID_old_to_new);
		for (int mm = V; mm < V + G; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_generated_group_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(generated_group_graph, vertexID_old_to_new);

		/*graph_v_of_v_idealID_DPBF_only_ec*/
		if (1) {
			double RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_v_of_v_idealID_DPBF_only_ec < < < 1 , 1 > > >(v_instance_graph, v_generated_group_graph, generated_group_vertices, RAM);
			cudaDeviceSynchronize()
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_DPBF_avg += (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_DPBF_sum += cost;

			if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}


		/*graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted*/
		if (1) {
			double RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_v_of_v_idealID_PrunedDPPlusPlus(v_instance_graph, v_generated_group_graph, generated_group_vertices, 1, RAM);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_PrunedDPPlusPlus_avg = time_PrunedDPPlusPlus_avg + (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_PrunedDPPlusPlus_sum = solution_cost_PrunedDPPlusPlus_sum + cost;

			if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}

		if (fabs(solution_cost_DPBF_sum-solution_cost_PrunedDPPlusPlus_sum) > 1e-8 ) {
			cout << "solution_cost_DPBF_sum=" << solution_cost_DPBF_sum << endl;
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