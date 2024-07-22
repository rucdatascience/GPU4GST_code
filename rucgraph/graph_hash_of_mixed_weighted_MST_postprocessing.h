
#pragma once


#include "graph_hash_of_mixed_weighted_minimum_spanning_tree.h"
#include "graph_hash_of_mixed_weighted_extract_subgraph_for_a_hash_of_vertices.h"



graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(graph_hash_of_mixed_weighted& sub_graph) {

	/*this function just returns the MST of the input sub_graph, the edges not in sub_graph is not considered;
	time complexity O(|subgraph_E|+|subgraph_V|log|subgraph_V|);
	nw and ec will be the same with sub_graph;
	*/

	graph_hash_of_mixed_weighted postprocessed_MST;

	for (auto it = sub_graph.hash_of_vectors.begin(); it != sub_graph.hash_of_vectors.end(); it++) {
		graph_hash_of_mixed_weighted_add_vertex(postprocessed_MST, it->first,
			sub_graph.hash_of_vectors[it->first].vertex_weight); // insert vertex
	}

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_hash_of_mixed_weighted_minimum_spanning_tree(sub_graph);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) { // predesscor of start_v is start_v
			graph_hash_of_mixed_weighted_add_edge(postprocessed_MST, v1, v2, graph_hash_of_mixed_weighted_edge_weight(sub_graph, v1, v2));
		}
	}

	return postprocessed_MST;

}



graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_MST_postprocessing(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& theta) {

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|+|subgraph_E|+|subgraph_V|log|subgraph_V|); nw and ec will be the same with input_graph, may not theta*/

	graph_hash_of_mixed_weighted postprocessed_theta; // 

	/*time complexity O(|subgraph_V|)*/
	std::unordered_set<int> hash_of_v;
	for (auto it = theta.hash_of_vectors.begin(); it != theta.hash_of_vectors.end(); it++) {
		hash_of_v.insert(it->first);
		graph_hash_of_mixed_weighted_add_vertex(postprocessed_theta, it->first,
			input_graph.hash_of_vectors[it->first].vertex_weight); // insert vertex
	}

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|)*/
	graph_hash_of_mixed_weighted subgraph_theta = graph_hash_of_mixed_weighted_extract_subgraph_for_a_hash_of_vertices(
		input_graph, hash_of_v);

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_hash_of_mixed_weighted_minimum_spanning_tree(subgraph_theta);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) { // predesscor of start_v is start_v
			graph_hash_of_mixed_weighted_add_edge
			(postprocessed_theta, v1, v2, graph_hash_of_mixed_weighted_edge_weight(input_graph, v1, v2));
		}
	}

	return postprocessed_theta;

}



graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_MST_postprocessing_hash_of_vertices
(graph_hash_of_mixed_weighted& input_graph, std::unordered_set<int>& hash_of_v) {

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|+|subgraph_E|+|subgraph_V|log|subgraph_V|); nw and ec in theta will be the same with input_graph*/

	graph_hash_of_mixed_weighted postprocessed_theta; // nw and ec will be the same with input_graph, may not theta

	/*time complexity O(|subgraph_V)*/
	for (auto it = hash_of_v.begin(); it != hash_of_v.end(); it++) {
		graph_hash_of_mixed_weighted_add_vertex(postprocessed_theta, *it, input_graph.hash_of_vectors[*it].vertex_weight); // insert vertex
	}

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|)*/
	graph_hash_of_mixed_weighted subgraph_theta = graph_hash_of_mixed_weighted_extract_subgraph_for_a_hash_of_vertices(
		input_graph, hash_of_v);

	//graph_hash_of_mixed_weighted_print(subgraph_theta);

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_hash_of_mixed_weighted_minimum_spanning_tree(subgraph_theta);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) {
			graph_hash_of_mixed_weighted_add_edge(postprocessed_theta, v1, v2, graph_hash_of_mixed_weighted_edge_weight(input_graph, v1, v2));
		}
	}

	return postprocessed_theta;

}






