#pragma once
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>

#include <vector>
#include "graph_v_of_v_idealID.h"
#include <unordered_map>
/*for GPU*/
using namespace std;

class CSR_graph
{
public:
    std::vector<int> INs_Neighbor_start_pointers, OUTs_Neighbor_start_pointers, ALL_start_pointers; // Neighbor_start_pointers[i] is the start point of neighbor information of vertex i in Edges and Edge_weights
    /*
        Now, Neighbor_sizes[i] = Neighbor_start_pointers[i + 1] - Neighbor_start_pointers[i].
        And Neighbor_start_pointers[V] = Edges.size() = Edge_weights.size() = the total number of edges.
    */
    std::vector<int> INs_Edges, OUTs_Edges, all_Edges;                      // Edges[Neighbor_start_pointers[i]] is the start of Neighbor_sizes[i] neighbor IDs
    std::vector<int> INs_Edge_weights, OUTs_Edge_weights, ALL_Edge_weights; // Edge_weights[Neighbor_start_pointers[i]] is the start of Neighbor_sizes[i] edge weights
    int *in_pointer, *out_pointer, *in_edge, *out_edge, *all_pointer, *all_edge;
    int *in_edge_weight, *out_edge_weight, *all_edge_weight;
    int E_all = 0, V;
    std::vector<int> pto;
    std::unordered_map<int, int> otn;
};

// CSR_graph<weight_type> toCSR(graph_structure<weight_type>& graph)
inline CSR_graph toCSR(graph_v_of_v_idealID graph)
{
//cudaSetDevice(1);
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
    cudaMallocManaged(&ARRAY.all_edge_weight, E_all * sizeof(int));

    cudaMemcpy(ARRAY.all_pointer, ARRAY.ALL_start_pointers.data(), (V + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ARRAY.all_edge, ARRAY.all_Edges.data(), E_all * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ARRAY.all_edge_weight, ARRAY.ALL_Edge_weights.data(), E_all * sizeof(int), cudaMemcpyHostToDevice);
     
    return ARRAY;
}
inline void toCSR_three(CSR_graph ARRAY[3], graph_v_of_v_idealID graph, int *community, int *c_size)
{

    for (size_t i = 0; i < 3; i++)
    {
        int V = c_size[i];
        ARRAY[i].V = V;
        ARRAY[i].INs_Neighbor_start_pointers.resize(V + 1); // Neighbor_start_pointers[V] = Edges.size() = Edge_weights.size() = the total number of edges.
        ARRAY[i].OUTs_Neighbor_start_pointers.resize(V + 1);
        ARRAY[i].ALL_start_pointers.resize(V + 1);
        ARRAY[i].pto.resize(V + 1);
    }
    int V = graph.size();
    int ptr[3] = {0, 0, 0}, acc[3] = {0, 0, 0}; // ptr pointer of vertex ; acc accord number of edges
    for (size_t i = 0; i < V; i++)
    {
        int c = community[i];
        ARRAY[c].pto[ptr[c]] = i; // to original vertex
        if (ARRAY[c].otn.find(i) == ARRAY[c].otn.end())
        {
            ARRAY[c].otn[i] = ptr[c]; // otn从原来指向现在坐标 pto指向原来坐标
        }
        ptr[c]++;
    }

    for (int i = 0; i < V; i++)
    {
        int c = community[i];
        ARRAY[c].ALL_start_pointers[ARRAY[c].otn[i]] = acc[c];
        for (int j = 0; j < graph[i].size(); j++)
        {
            if (community[graph[i][j].first] == c)
            {
                ARRAY[c].all_Edges.push_back(ARRAY[c].otn[graph[i][j].first]);
                ARRAY[c].ALL_Edge_weights.push_back(graph[i][j].second);
                ARRAY[c].E_all++;
                acc[c]++;
            }
        }
    }
    std::cout << "111" << std::endl;
    for (size_t i = 0; i < 3; i++)
    {
        ARRAY[i].ALL_start_pointers[ARRAY[i].V] = acc[i];
        int E_all = acc[i];
        V = ARRAY[i].V;

        cudaMallocManaged(&ARRAY[i].all_pointer, (V + 1) * sizeof(int));
        cudaMallocManaged(&ARRAY[i].all_edge, E_all * sizeof(int));
        cudaMallocManaged(&ARRAY[i].all_edge_weight, E_all * sizeof(int));

        cudaMemcpy(ARRAY[i].all_pointer, ARRAY[i].ALL_start_pointers.data(), (V + 1) * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(ARRAY[i].all_edge, ARRAY[i].all_Edges.data(), E_all * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(ARRAY[i].all_edge_weight, ARRAY[i].ALL_Edge_weights.data(), E_all * sizeof(int), cudaMemcpyHostToDevice);
    }
    /*for (size_t k = 0; k < 3; k++)
     {cout << "graph part" << k << endl;
         for (size_t i = 0; i < ARRAY[k].V; i++)
         {

             for (size_t j = ARRAY[k].all_pointer[i]; j < ARRAY[k].all_pointer[i + 1]; j++)
             {
                 cout << "edge to " << ARRAY[k].all_edge[j] << " weight " << ARRAY[k].all_edge_weight[j] << " ";
             }
             cout << endl;
         }
     }*/
}