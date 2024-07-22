#include <chrono>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.hpp>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_connected_components.hpp>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_connected_graph.hpp>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.hpp>
#include <graph_hash_of_mixed_weighted_generate_random_groups_of_vertices.hpp>
#include <graph_hash_of_mixed_weighted_save_for_GSTP.hpp>
#include <graph_hash_of_mixed_weighted_read_for_GSTP.hpp>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID.hpp>
#include <graph_hash_of_mixed_weighted_sum_of_nw_ec.hpp>
#include <graph_v_of_v_idealID_DPBF_only_ec.cu>
#include <graph_v_of_v_idealID_PrunedDPPlusPlus.hpp>

bool this_is_a_feasible_solution(graph_hash_of_mixed_weighted& solu, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& group_vertices);

extern "C" 

void test_graph_v_of_v_idealID_DPBF_only_ec();