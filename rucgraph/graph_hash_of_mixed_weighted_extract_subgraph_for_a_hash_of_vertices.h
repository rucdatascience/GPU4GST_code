#pragma once


graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_extract_subgraph_for_a_hash_of_vertices(
	graph_hash_of_mixed_weighted& input_graph, std::unordered_set<int>& hash_of_v) {

	/*extract a smaller_graph, which contains all the vertices in hash_of_v,
	and all the edges between vertices in hash_of_v;
	time complexity O(|V_list|+|adj_v of V_list in input_graph|)*/

	/*graph_unordered_map_extract_subgraph_for_a_list_of_vertices.h first transform a list to a hash,
	and then do the following*/


	/*time complexity O(|V_list|+|adj_v of V_list in input_graph|)*/
	graph_hash_of_mixed_weighted smaller_graph;
	for (auto i = hash_of_v.begin(); i != hash_of_v.end(); i++) {
		int v1 = *i;
		int nw = input_graph.hash_of_vectors[v1].vertex_weight;
		graph_hash_of_mixed_weighted_add_vertex(smaller_graph, v1, nw); // add vertex

		auto search = input_graph.hash_of_hashs.find(v1);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int v2 = it2->first;
				int ec = it2->second;
				if (hash_of_v.count(v2) > 0 & v2 > v1) { // v2 is in the list and only add edge once
					graph_hash_of_mixed_weighted_add_edge(smaller_graph, v1, v2, ec); // add edge
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(v1);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int v2 = it2->first;
				int ec = it2->second;
				if (hash_of_v.count(v2) > 0 & v2 > v1) { // v2 is in the list and only add edge once
					graph_hash_of_mixed_weighted_add_edge(smaller_graph, v1, v2, ec); // add edge
				}
			}
		}
	}

	return smaller_graph;

}
