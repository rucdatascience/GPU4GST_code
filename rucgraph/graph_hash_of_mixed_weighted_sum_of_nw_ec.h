#pragma once


int graph_hash_of_mixed_weighted_sum_of_nw_ec(graph_hash_of_mixed_weighted& input_graph) {

	int sum = 0;

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		int w_i = it1->second.vertex_weight;
		sum = sum + w_i;

		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					int c_ij = it2->second;
					sum = sum + c_ij;
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					int c_ij = it2->second;
					sum = sum + c_ij;
				}
			}
		}
	}

	return sum;

}




int graph_hash_of_mixed_weighted_sum_of_nw(graph_hash_of_mixed_weighted& input_graph) {

	int sum = 0;

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		int w_i = it1->second.vertex_weight;
		sum = sum + w_i;
	}

	return sum;

}



int graph_hash_of_mixed_weighted_sum_of_ec(graph_hash_of_mixed_weighted& input_graph) {

	int sum = 0;

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;

		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					sum = sum + it2->second;
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					sum = sum + it2->second;
				}
			}
		}
	}

	return sum;

}

