
#include <DP.cuh>
using namespace std;
typedef struct node
{
    int type;    // =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge
    double cost; // cost of this tree T(v,p);
    int u;       // if this tree is built by grown, then it's built by growing edge (v,u);
    int p1, p2;  // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);
} node;
int E, N, problem_size, width, height;
int *all_pointer, *all_edge, *non_overlapped_group_sets_IDs_pointer, *non_overlapped_group_sets_IDs_gpu, *tes;
dim3 blockPerGrid, threadPerGrid;
double *edge_cost, inf = 10240;
node *tree;
int graph_v_of_v_idealID_DPBF_vertex_group_set_ID_gpu(int vertex, graph_v_of_v_idealID &group_graph,
                                                      std::unordered_set<int> &cumpulsory_group_vertices)
{

    /*time complexity: O(|Gamma|); this function returns the maximum group set ID for a single vertex*/
    // if group i have edge to v,v will give bit i value 1;
    int ID = 0;
    int pow_num = 0;
    for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++)
    {
        if (graph_v_of_v_idealID_contain_edge(group_graph, vertex, *it))
        { // vertex is in group *it
            ID = ID + pow(2, pow_num);
        }
        pow_num++;
    }

    return ID;
}
__global__ void test_and(int *tes, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        tes[idx] = 6 ^ idx;
    }
}
__global__ void init(node *tree, int range, size_t pitch_node, double inf, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        node *row = (node *)((char *)tree + idx * pitch_node);
        for (int subs = 1; subs <= range; subs++)
            row[subs].cost = inf;
    }
}
__global__ void two(node *tree, size_t pitch_node, int s, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        node *row = (node *)((char *)tree + idx * pitch_node);
        for (size_t i = 0; i < 5; i++)
        {
            row[i].cost = i;
        }
    }
}
__global__ void merge_gpu(node *tree, size_t pitch_node, int s, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        node *row = (node *)((char *)tree + idx * pitch_node);
        for (int subs = s & (s - 1); subs; subs = s & (subs - 1))
            if (row[s].cost > row[subs].cost + row[s ^ subs].cost)
            {
                row[s].cost = row[subs].cost + row[s ^ subs].cost;
                row[s].type = 2;
                row[s].p1 = subs;
                row[s].p2 = s ^ subs;
            }
    }
}

__global__ void bell_ford(node *tree, int *edge, int *edge_pointer, size_t pitch_node, double *edge_cost, int s, int N)
{
    int v = blockIdx.x * blockDim.x + threadIdx.x;
    if (v < N)
    {
        node *row_v = (node *)((char *)tree + v * pitch_node);
        for (int i = edge_pointer[v]; i < edge_pointer[v + 1]; i++)
        {
            int u = edge[i];
            node *row_u = (node *)((char *)tree + u * pitch_node);

            if (row_v[s].cost > row_u[s].cost + edge_cost[i])
            {
                row_v[s].cost = row_u[s].cost + edge_cost[i];
                row_v[s].type = 1;
                row_v[s].u = u;
            }
        }
    }
}
graph_hash_of_mixed_weighted DPBF(CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph)
{
    E = graph.E_all;
    N = graph.V;
    all_edge = graph.all_edge, all_pointer = graph.all_pointer, edge_cost = graph.all_edge_weight;

    int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1;
    cudaMallocManaged((void **)&tes, sizeof(int) * 10);
    test_and<<<1, 10>>>(tes, 10);
    cudaDeviceSynchronize();

    cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_pointer, sizeof(int) * (group_sets_ID_range + 2));
    std::cout << "group range " << group_sets_ID_range << std::endl;
    threadPerGrid.x = THREAD_PER_BLOCK, blockPerGrid.x = (N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK;
    width = group_sets_ID_range + 1, height = N;
    size_t pitch_node;
    node host_tree[height][width];
    // init<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, group_sets_ID_range, pitch_node, inf,N);
    cudaMallocPitch(&tree, &pitch_node, width * sizeof(node), height);
    std::cout << "pitch " << pitch_node << " width " << width << std::endl;
    for (int v = 0; v < N; v++)
    {
        int group_set_ID_v = graph_v_of_v_idealID_DPBF_vertex_group_set_ID_gpu(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
        // cout<<group_set_ID_v<<" for "<<v<<endl;
        for (int p = 1; p <= group_sets_ID_range; p++)

        { // p is non-empty; time complexity: O(2^|Gamma|) //get all its subset ,which is required in next merge and grow steps
            host_tree[v][p].cost = inf;
            if ((p | group_set_ID_v) == group_set_ID_v)
            { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v
                /*T(v,p)*/
                node node;
                node.cost = 0;
                node.type = 0;
                host_tree[v][p] = node;
            }
        }
    }
    std::cout << "cost init " << std::endl;
    cudaMemcpy2D(tree, pitch_node, host_tree, width * sizeof(node), width * sizeof(node), height, cudaMemcpyHostToDevice);
    cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);

    cout << "after init " << endl;

    cudaDeviceSynchronize();
    // Device code
    for (int r = 1; r <= group_sets_ID_range; r++)
    {

        merge_gpu<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, pitch_node, r, N);
        cudaDeviceSynchronize();
        for (int round = 1; round < N; round++)
        {
            bell_ford<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, all_edge, all_pointer, pitch_node, edge_cost, r, N);
        cudaDeviceSynchronize();
        }

        // two<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, pitch_node, r, N);
        
    }

    std::cout << "while over" << std::endl;
    cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);

    std::cout << "all copy complete" << std::endl;

    int min_cost = 1e9, min_node = -1;
    for (int i = 0; i < N; i++)
    {
        if (host_tree[i][group_sets_ID_range].cost < min_cost)
        {
            min_cost = host_tree[i][group_sets_ID_range].cost;
            min_node = i;
        }
    }

    std::cout << "root " << min_node << std::endl;
    graph_hash_of_mixed_weighted solution_tree;
    std::queue<std::pair<int, int>> waited_to_processed_trees; // <v, p>
    int root_v = min_node, root_p = group_sets_ID_range;
    waited_to_processed_trees.push({root_v, root_p});

    while (waited_to_processed_trees.size() > 0)
    {

        int v = waited_to_processed_trees.front().first, p = waited_to_processed_trees.front().second;
        waited_to_processed_trees.pop();

        graph_hash_of_mixed_weighted_add_vertex(solution_tree, v, 0);

        auto pointer_trees_v_p = host_tree[v][p];
        int form_type = pointer_trees_v_p.type;
        if (form_type == 0)
        { // T(v,p) is a single vertex
        }
        else if (form_type == 1)
        { // T(v,p) is formed by grow
            int u = host_tree[v][p].u;

            waited_to_processed_trees.push({u, p});
            /*insert (u,v); no need to insert weight of u here, which will be inserted later for T(u,p)*/
            double c_uv = graph_v_of_v_idealID_edge_weight(input_graph, u, v);
            graph_hash_of_mixed_weighted_add_edge(solution_tree, u, v, c_uv);
        }
        else
        { // T(v,p) is formed by merge
            int p1 = host_tree[v][p].p1, p2 = host_tree[v][p].p2;

            waited_to_processed_trees.push({v, p1});
            waited_to_processed_trees.push({v, p2});
        }
    }
    return solution_tree;
}
