#pragma once
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>

#include <vector>
#include "graph_v_of_v_idealID.h"

/*for GPU*/

class CSR_graph
{
public:
    std::vector<int> INs_Neighbor_start_pointers, OUTs_Neighbor_start_pointers, ALL_start_pointers; // Neighbor_start_pointers[i] is the start point of neighbor information of vertex i in Edges and Edge_weights
    /*
        Now, Neighbor_sizes[i] = Neighbor_start_pointers[i + 1] - Neighbor_start_pointers[i].
        And Neighbor_start_pointers[V] = Edges.size() = Edge_weights.size() = the total number of edges.
    */
    std::vector<int> INs_Edges, OUTs_Edges, all_Edges;                              // Edges[Neighbor_start_pointers[i]] is the start of Neighbor_sizes[i] neighbor IDs
    std::vector<double> INs_Edge_weights, OUTs_Edge_weights, ALL_Edge_weights; // Edge_weights[Neighbor_start_pointers[i]] is the start of Neighbor_sizes[i] edge weights
    int *in_pointer, *out_pointer, *in_edge, *out_edge, *all_pointer, *all_edge;
    double *in_edge_weight, *out_edge_weight, *all_edge_weight;
    int E_all=0, V;
};


// CSR_graph<weight_type> toCSR(graph_structure<weight_type>& graph)
inline CSR_graph toCSR(graph_v_of_v_idealID graph)
{

    CSR_graph ARRAY;

    int V = graph.size();
    ARRAY.V = V;
    ARRAY.INs_Neighbor_start_pointers.resize(V + 1); // Neighbor_start_pointers[V] = Edges.size() = Edge_weights.size() = the total number of edges.
    ARRAY.OUTs_Neighbor_start_pointers.resize(V + 1);
    ARRAY.ALL_start_pointers.resize(V + 1);

    int pointer = 0;
    for (int i = 0; i < V; i++)
    {
        ARRAY.ALL_start_pointers[i] = pointer;
        for (int j = 0; j < graph[i].size(); j++)
        {
            ARRAY.all_Edges.push_back(graph[i][j].first);
            ARRAY.ALL_Edge_weights.push_back(graph[i][j].second);
            ARRAY.E_all++;
        }

        pointer += graph[i].size();
    }
    ARRAY.ALL_start_pointers[V] = pointer;

    int E_all = ARRAY.E_all;

    cudaMallocManaged(&ARRAY.all_pointer, (V + 1) * sizeof(int));
    cudaMallocManaged(&ARRAY.all_edge, E_all * sizeof(int));
    cudaMallocManaged(&ARRAY.all_edge_weight, E_all * sizeof(double));

    cudaMemcpy(ARRAY.all_pointer, ARRAY.ALL_start_pointers.data(), (V + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ARRAY.all_edge, ARRAY.all_Edges.data(), E_all * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ARRAY.all_edge_weight, ARRAY.ALL_Edge_weights.data(), E_all * sizeof(double), cudaMemcpyHostToDevice);
    return ARRAY;
}
