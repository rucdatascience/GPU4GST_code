#ifndef DP
#define DP
#pragma once
#include <iostream>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <stdio.h>
#include <string>
#include <ctime>
#include <graph_v_of_v_idealID/csr_graph.hpp>
#include "../rucgraph/graph_v_of_v_idealID/csr_graph.hpp"
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>


#define THREAD_PER_BLOCK 512
graph_hash_of_mixed_weighted DPBF(CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph);

#endif