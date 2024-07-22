#include <chrono>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include "graph_hash_of_mixed_weighted_generate_random_groups_of_vertices.hpp"
#include "graph_hash_of_mixed_weighted_save_for_GSTP.hpp"
#include "graph_hash_of_mixed_weighted_read_for_GSTP.hpp"
#include "graph_v_of_v_idealID/graph_v_of_v_idealID.hpp"
#include "graph_hash_of_mixed_weighted_sum_of_nw_ec.hpp"

void getpath(int statu,int x,graph_hash_of_mixed_weighted& solu);

extern "C" 

__global__ graph_hash_of_mixed_weighted graph_v_of_v_idealID_DPBF_only_ec(
	graph_v_of_v_idealID& input_graph, graph_v_of_v_idealID& group_graph, std::unordered_set<int>& cumpulsory_group_vertices, double& RAM_MB);