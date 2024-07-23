
#include <DPQ.cuh>
#include <cuco/static_set.cuh>
#include <thrust/device_vector.h>
using namespace std;
/*this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.

time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )*/

struct DPBF_min_node
{
	int v;
	int p; // group_set_ID
	float priority_value;
};
typedef struct node
{
	int type;	// =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge
	int cost;	// cost of this tree T(v,p);
	int u;		// if this tree is built by grown, then it's built by growing edge (v,u);
	int p1, p2; // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);
} node;
int E, N, width, height;
int *type, *u, *p1, *p2, *visit, *queue_size, *tes, *queue_end;
int *all_pointer, *all_edge, *edge_cost, *non_overlapped_group_sets_IDs_pointer, *non_overlapped_group_sets_IDs_gpu;
dim3 blockPerGrid, threadPerGrid;
node *tree;
int *Queue;
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
template <typename T>
graph_hash_of_mixed_weighted graph_v_of_v_idealID_DPBF_build_tree_gpu(int root_v, int root_p, graph_v_of_v_idealID &input_graph,
																	  T trees)
{

	/*this function builds tree T(v,p) at a cost of O(|V|)*/

	graph_hash_of_mixed_weighted solution_tree;

	std::queue<std::pair<int, int>> waited_to_processed_trees; // <v, p>
	waited_to_processed_trees.push({root_v, root_p});

	while (waited_to_processed_trees.size() > 0)
	{

		int v = waited_to_processed_trees.front().first, p = waited_to_processed_trees.front().second;
		waited_to_processed_trees.pop();

		/*insert v*/
		graph_hash_of_mixed_weighted_add_vertex(solution_tree, v, 0);

		auto pointer_trees_v_p = trees[v][p];
		int form_type = pointer_trees_v_p.type;
		if (form_type == 0)
		{ // T(v,p) is a single vertex
		}
		else if (form_type == 1)
		{ // T(v,p) is formed by grow
			int u = trees[v][p].u;
			waited_to_processed_trees.push({u, p});
			/*insert (u,v); no need to insert weight of u here, which will be inserted later for T(u,p)*/
			int c_uv = graph_v_of_v_idealID_edge_weight(input_graph, u, v);
			graph_hash_of_mixed_weighted_add_edge(solution_tree, u, v, c_uv);
		}
		else
		{ // T(v,p) is formed by merge
			int p1 = trees[v][p].p1, p2 = trees[v][p].p2;
			waited_to_processed_trees.push({v, p1});
			waited_to_processed_trees.push({v, p2});
		}
	}

	return solution_tree;
}
void graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(int group_sets_ID_range)
{

	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID;

	time complexity: O(4^|Gamma|), since group_sets_ID_range=2^|Gamma|;

	the original DPBF code use the same method in this function, and thus has the same O(4^|Gamma|) complexity;*/

	std::vector<int> non_overlapped_group_sets_IDs; // <set_ID, non_overlapped_group_sets_IDs>
	int len = 0;
	for (int i = 1; i <= group_sets_ID_range; i++)
	{ // i is a nonempty group_set ID
		non_overlapped_group_sets_IDs_pointer[i] = len;
		for (int j = 1; j < group_sets_ID_range; j++)
		{ // j is another nonempty group_set ID
			if ((i & j) == 0)
			{ // i and j are non-overlapping group sets
				/* The & (bitwise AND) in C or C++ takes two numbers as operands and does AND on every bit of two numbers. The result of AND for each bit is 1 only if both bits are 1.
				https://www.programiz.com/cpp-programming/bitwise-operators */
				non_overlapped_group_sets_IDs.push_back(j);
				len++;
			}
		}
	}
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_gpu, sizeof(int) * len);
	cudaMemcpy(non_overlapped_group_sets_IDs_gpu, non_overlapped_group_sets_IDs.data(), sizeof(int) * len, cudaMemcpyHostToDevice);
	std::cout << "len= " << len << std::endl;
}

__global__ void reset(int *visit, int N, size_t pitch_vis, int problem_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		int *row = (int *)((char *)visit + idx * pitch_vis);
		for (int i = 0; i <= problem_size; i++)
		{
			row[i] = -1;
		}
	}
}
template <typename SetRef>
__global__ void Relax(int *Queue, int *queue_size, int *sets_IDs, int *sets_IDS_pointer, int *edge, int *edge_cost, int *pointer, size_t pitch_node, node *tree, int width, SetRef set)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *queue_size)
	{
		int top_node = Queue[idx];
		int v = top_node / width, p = top_node % width;
		namespace cg = cooperative_groups;

		constexpr auto cg_size = SetRef::cg_size;

		auto tile = cg::tiled_partition<cg_size>(cg::this_thread_block());
		/*grow*/
		for (int i = pointer[v]; i < pointer[v + 1]; i++)
		{

			int u = edge[i];
			int cost_euv = edge_cost[i];
			node *row = (node *)((char *)tree + v * pitch_node);
			int grow_tree_cost = cost_euv + row[p].cost;
			row = (node *)((char *)tree + u * pitch_node);
			int old = atomicMin(&row[p].cost, grow_tree_cost);
			if (old <= grow_tree_cost)
				continue;

			else
			{
				/*below: O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since each u is checked 2^|Gamma| times, and Q_T contains at most 2^|Gamma|*V elements */

				/*update T(u,p) by grow T(v,p) with (u,v)*/
				row[p].type = 1;
				row[p].u = v;
				set.insert(tile, u * width + p);

				/* int *vis_row = (int *)((char *)visit + u * pitch_vis);
				if (vis_row[p] == -1)
				{ // T(u,p) is not in Q_T
					int pos = atomicAdd(next_queue_size, 1);
					next_queue[pos] = x;
					vis_row[p] = pos;
				}
				else
				{ // T(u,p) is in Q_T
					next_queue[vis_row[p]] = x;
				} */
			}
		}

		/*merge*/
		int p1 = p;
		for (auto it = sets_IDS_pointer[p1]; it < sets_IDS_pointer[p1 + 1]; it++)
		{
			int p2 = sets_IDs[it]; // p2 is not overlapped with p1
			int cost_Tvp1, cost_Tvp2;
			node *row = (node *)((char *)tree + v * pitch_node);
			cost_Tvp1 = row[p1].cost;
			cost_Tvp2 = row[p2].cost;
			int p1_cup_p2 = p1 + p2;

			int merged_tree_cost = cost_Tvp1 + cost_Tvp2;
			int old = atomicMin(&row[p1_cup_p2].cost, merged_tree_cost);
			if (old > merged_tree_cost)
			{ // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				row[p1_cup_p2].type = 2;
				row[p1_cup_p2].p1 = p1;
				row[p1_cup_p2].p2 = p2;
				set.insert(tile, v * width + p1_cup_p2);
			}
		}
	}
}

graph_hash_of_mixed_weighted DPBF_GPU(CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph)
{
	E = graph.E_all;
	N = graph.V;
	all_edge = graph.all_edge, all_pointer = graph.all_pointer, edge_cost = graph.all_edge_weight;
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1;
	long long unsigned int problem_size = N * pow(2, cumpulsory_group_vertices.size());
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_pointer, sizeof(int) * (group_sets_ID_range + 2));
	cudaMallocManaged((void **)&queue_size, sizeof(int));
	cudaMallocManaged((void **)&Queue, sizeof(int) * problem_size);
	cudaMallocManaged((void **)&tes, sizeof(int) * THREAD_PER_BLOCK);
	graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(group_sets_ID_range);
	std::cout << "group range " << group_sets_ID_range << std::endl;
	// thrust::device_vector<node> tree(problem_size);
	threadPerGrid.x = THREAD_PER_BLOCK;
	blockPerGrid.x = (N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK;
	width = group_sets_ID_range + 1, height = N;
	size_t pitch_node;
	node host_tree[height][width];
	cudaMallocPitch(&tree, &pitch_node, width * sizeof(node), height);
	std::cout << "pitch " << pitch_node << " " << " width " << width << std::endl;
	auto constexpr load_factor = 0.5;
	using Key = int;

	int constexpr empty_key_sentinel = -1;
	std::size_t const capacity = std::ceil(problem_size / load_factor);
	cuco::static_set<Key> set{problem_size, cuco::empty_key{empty_key_sentinel}};
	*queue_size = 0;
	for (int v = 0; v < N; v++)
	{	host_tree[v][0].cost = 0;
		int group_set_ID_v = graph_v_of_v_idealID_DPBF_vertex_group_set_ID_gpu(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_sets_ID_range; p++)
		{ // p is non-empty; time complexity: O(2^|Gamma|) //get all its subset ,which is required in next merge and grow steps
			host_tree[v][p].cost = 1024;

			if ((p | group_set_ID_v) == group_set_ID_v)
			{ // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v
				/*T(v,p)*/
				node node;
				node.cost = 0;
				node.type = 0;
				host_tree[v][p] = node;

				Queue[*queue_size] = v * width + p;
				*queue_size += 1;
				/* DPBF_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0; */
			}
		}
	}
	/* 	std::cout << "queue init " << std::endl;
		for (size_t i = 0; i < 100; i++)
		{
			std::cout << " v " << Queue[i].v << " p " << Queue[i].p << "; ";
		} */
	set.insert(Queue, Queue + *queue_size);
	queue_end = set.retrieve_all(Queue);
	cout<<"init set size"<<(queue_end - Queue)<<endl;
	cudaMemcpy2D(tree, pitch_node, host_tree, width * sizeof(node), width * sizeof(node), height, cudaMemcpyHostToDevice);
	std::cout << "queue size init " << *queue_size << std::endl;
	for (size_t i = 0; i < 10; i++)
	{
		for (size_t j = 0; j <= group_sets_ID_range; j++)
		{
			cout << host_tree[i][j].cost << " ";
		}
		cout << endl;
	}
	int r = 0;
	while (*queue_size != 0)
	{
		std::cout << "round " << r++ << std::endl;
		std::cout << "queue size " << *queue_size << std::endl;
		for (size_t i = 0; i < *queue_size; i++)
		{
			std::cout << " v " << Queue[i] / group_sets_ID_range << " p " << Queue[i] % group_sets_ID_range << "; ";
		}
		// reset<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(visit, N, pitch_vis, group_sets_ID_range);
		Relax<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(Queue, queue_size, non_overlapped_group_sets_IDs_gpu,
																							 non_overlapped_group_sets_IDs_pointer, all_edge, edge_cost, all_pointer, pitch_node, tree, width, set.ref(cuco::insert));
		cudaDeviceSynchronize();
		cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
/* 		for (size_t i = 0; i < 10; i++)
		{
			for (size_t j = 0; j <= group_sets_ID_range; j++)
			{
				cout << host_tree[i][j].cost << " ";
			}
			cout << endl;
		} */
		queue_end = set.retrieve_all(Queue);
		*queue_size = (queue_end - Queue);
		set.clear();
	}
	std::cout << "while over" << std::endl;
	cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);

	std::cout << "all copy complete ,now list cost " << std::endl;
	int min_cost = 1e9, min_node = -1;
	for (int i = 0; i < N; i++)
	{
		cout << host_tree[i][group_sets_ID_range].cost << " ";
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
			int c_uv = graph_v_of_v_idealID_edge_weight(input_graph, u, v);
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
