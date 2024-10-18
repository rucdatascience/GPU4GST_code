
#include <DPQ.cuh>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <omp.h>
using namespace std;
// use one label lowerbound
/*this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.
time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )*/

typedef struct queue_element
{
	int v, p, cost;
	bool operator<(const queue_element &other) const
	{
		return cost < other.cost;
	}
} queue_element;

std::vector<int> non_overlapped_group_sets_IDs_pointer_host, non_overlapped_group_sets_IDs;
void set_max_ID(graph_v_of_v_idealID &group_graph, std::unordered_set<int> &cumpulsory_group_vertices, node **host_tree)
{
	int bit_num = 1, v;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++, bit_num <<= 1)
	{

		for (size_t to = 0; to < group_graph[*it].size(); to++)
		{
			v = group_graph[*it][to].first;
			host_tree[v][bit_num].cost = 0;
		}
	}
}
int get_max(int vertex, node **host_tree, int G)
{
	int re = 0;
	for (size_t i = 1; i < G; i <<= 1)
	{
		if (host_tree[vertex][i].cost == 0)
		{
			re += i;
		}
	}

	return re;
}
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

__global__ void Relax_1(queue_element *Queue_dev, int *queue_size, queue_element *new_queue_device, int *new_queue_size, int *sets_IDs, int *sets_IDS_pointer, int *edge, int *edge_cost, int *pointer, int width, node *tree, int inf, int *best, int full, int *lb0)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *queue_size)
	{
		queue_element top_node = Queue_dev[idx];
		int v = top_node.v, p = top_node.p;
		int x_slash = full - p;
		tree[v * width + p].update = 0;
		if (tree[v * width + x_slash].cost != inf)
		{
			int new_best = tree[v * width + x_slash].cost + tree[v * width + p].cost;
			atomicMin(best, new_best);
			if (new_best < *best)
			{
				atomicMin(&tree[v * width + full].cost, new_best);
				tree[v * width + full].type = 2;
				tree[v * width + full].p1 = p;
				tree[v * width + full].p2 = x_slash;
			}
		}
		else
		{
			int new_best = tree[v * width + p].cost;
			for (size_t i = 1; i <= x_slash; i <<= 1)
			{
				if (i & x_slash)
				{
					new_best += tree[v * width + i].cost;
				}
			}
			atomicMin(best, new_best);
		}
		if (tree[v * width + p].cost > (*best) / 2)
		{
			return;
		}
		for (int i = pointer[v]; i < pointer[v + 1]; i++)
		{
			/*grow*/

			int u = edge[i];
			int cost_euv = edge_cost[i];
			int grow_tree_cost = tree[v * width + p].cost + cost_euv;
			int old = atomicMin(&tree[u * width + p].cost, grow_tree_cost);

			if (old >= grow_tree_cost && grow_tree_cost != inf)
			{
				tree[u * width + p].type = 1;
				tree[u * width + p].u = v;
				// enqueue operation
				int check = atomicCAS(&tree[u * width + p].update, 0, 1);
				if (check == 0 && lb0[u] + grow_tree_cost <= (*best))
				{
					int pos = atomicAdd(new_queue_size, 1);
					new_queue_device[pos].v = u;
					new_queue_device[pos].p = p;
					new_queue_device[pos].cost = grow_tree_cost;
				}
			}
		}

		/*merge*/
		int p1 = p;
		for (auto it = sets_IDS_pointer[p1]; it < sets_IDS_pointer[p1 + 1]; it++)
		{
			int p2 = sets_IDs[it]; // p2 is not overlapped with p1
			int cost_Tvp1 = tree[v * width + p1].cost, cost_Tvp2 = tree[v * width + p2].cost;
			int p1_cup_p2 = p1 + p2;
			int merged_tree_cost = cost_Tvp1 + cost_Tvp2;
			// && merged_tree_cost < 2 / 3 * (*best)
			int old = atomicMin(&tree[v * width + p1_cup_p2].cost, merged_tree_cost);

			if (old >= merged_tree_cost && merged_tree_cost != inf)
			{ // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				tree[v * width + p1_cup_p2].type = 2;
				tree[v * width + p1_cup_p2].p1 = p1;
				tree[v * width + p1_cup_p2].p2 = p2;

				if (merged_tree_cost < 0.667 * (*best) && lb0[v] + merged_tree_cost < (*best))
				{
					int check = atomicCAS(&tree[v * width + p1_cup_p2].update, 0, 1);
					if (check == 0)
					{
						int pos = atomicAdd(new_queue_size, 1);
						new_queue_device[pos].v = v;
						new_queue_device[pos].p = p1_cup_p2;
						new_queue_device[pos].cost = merged_tree_cost;
					}
				}
			}
		}
	}
}
__global__ void pre_scan(queue_element *new_queue_device, int *new_queue_size, queue_element *queue_device, int *queue_size, int *best, int up, double percent, int *check, queue_element *mid_queue, int *mid_queue_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *new_queue_size)
	{
		int cost = new_queue_device[idx].cost;
		if (cost < *best * percent)
		{
			int pos = atomicAdd(queue_size, 1);
			if (pos < up)
			{
				check[idx] = 1;
				queue_device[pos] = new_queue_device[idx];
			}
		}
		if (cost > *best)
		{
			check[idx] = 1;
		}
		if (check[idx] == 0)
		{
			int pos = atomicAdd(mid_queue_size, 1);
			mid_queue[pos] = new_queue_device[idx];
		}
	}
}

__global__ void BFS(int *bfs_queue, int *bfs_size, int *new_bfs_queue, int *new_bfs_size, int *edge, int *visited, int *pointer)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *bfs_size)
	{
		int v = bfs_queue[idx], pos;
		visited[v] = 1;
		for (size_t i = pointer[v]; i < pointer[v + 1]; i++)
		{
			int u = edge[i];
			if (visited[u] == 0)
			{
				pos = atomicAdd(new_bfs_size, 1);
				new_bfs_queue[pos] = u;
				visited[u] = 1;
			}
		}
	}
}

__global__ void Relax(queue_element *Queue_dev, int *queue_size, queue_element *new_queue_device, int *new_queue_size, int *sets_IDs, int *sets_IDS_pointer, int *edge, int *edge_cost, int *pointer, size_t pitch_node, node *tree, int inf, int *best, int full, int *lb0)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *queue_size)
	{
		queue_element top_node = Queue_dev[idx];
		int v = top_node.v, p = top_node.p;
		node *row_node_v = (node *)((char *)tree + v * pitch_node);
		int x_slash = full - p;
		row_node_v[p].update = 0;
		if (row_node_v[x_slash].cost != inf)
		{
			int new_best = row_node_v[x_slash].cost + row_node_v[p].cost;
			atomicMin(best, new_best);
			if (new_best < *best)
			{
				atomicMin(&row_node_v[full].cost, new_best);
				row_node_v[full].type = 2;
				row_node_v[full].p1 = p;
				row_node_v[full].p2 = x_slash;
			}
		}
		else
		{
			int new_best = row_node_v[p].cost;
			for (size_t i = 1; i <= x_slash; i <<= 1)
			{
				if (i & x_slash)
				{
					new_best += row_node_v[i].cost;
				}
			}
			atomicMin(best, new_best);
		}
		if (row_node_v[p].cost >= (*best) / 2)
		{
			return;
		}
		for (int i = pointer[v]; i < pointer[v + 1]; i++)
		{
			/*grow*/

			int u = edge[i];
			int cost_euv = edge_cost[i];
			node *row_u = (node *)((char *)tree + u * pitch_node);
			int grow_tree_cost = row_node_v[p].cost + cost_euv;
			int old = atomicMin(&row_u[p].cost, grow_tree_cost);
			node *row_node_u = (node *)((char *)tree + u * pitch_node);
			if (old >= grow_tree_cost && grow_tree_cost != inf)
			{
				row_node_u[p].type = 1;
				row_node_u[p].u = v;
				// enqueue operation
				int check = atomicCAS(&row_node_u[p].update, 0, 1);
				if (check == 0 && lb0[u] + grow_tree_cost <= (*best))
				{
					int pos = atomicAdd(new_queue_size, 1);
					new_queue_device[pos].v = u;
					new_queue_device[pos].p = p;
				}
			}
		}

		/*merge*/
		int p1 = p;
		for (auto it = sets_IDS_pointer[p1]; it < sets_IDS_pointer[p1 + 1]; it++)
		{
			int p2 = sets_IDs[it]; // p2 is not overlapped with p1
			int cost_Tvp1 = row_node_v[p1].cost, cost_Tvp2 = row_node_v[p2].cost;
			int p1_cup_p2 = p1 + p2;
			int merged_tree_cost = cost_Tvp1 + cost_Tvp2;
			// && merged_tree_cost < 2 / 3 * (*best)
			int old = atomicMin(&row_node_v[p1_cup_p2].cost, merged_tree_cost);

			if (old >= merged_tree_cost && merged_tree_cost != inf)
			{ // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				row_node_v[p1_cup_p2].type = 2;
				row_node_v[p1_cup_p2].p1 = p1;
				row_node_v[p1_cup_p2].p2 = p2;

				if (merged_tree_cost < 0.667 * (*best) && lb0[v] + merged_tree_cost < (*best))
				{
					int check = atomicCAS(&row_node_v[p1_cup_p2].update, 0, 1);
					if (check == 0)
					{
						int pos = atomicAdd(new_queue_size, 1);
						new_queue_device[pos].v = v;
						new_queue_device[pos].p = p1_cup_p2;
					}
				}
			}
		}
	}
}
__global__ void one_label_lb_1(node *tree, int width, int *lb0, int N, int G, int inf, int *can_find_solution)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N)
	{
		int row = idx * width;

		for (size_t i = 1; i < G; i <<= 1)
		{
			if (lb0[idx] < tree[row + i].cost)
			{
				lb0[idx] = tree[row + i].cost;
			}
		}
		if (lb0[idx] != inf)
		{
			*can_find_solution = 1;
		}

		for (size_t i = 1; i < G; i <<= 1)
		{
			for (size_t j = i << 1; j < G; j <<= 1)
			{
				if (tree[row + i + j].cost < tree[row + i].cost + tree[row + j].cost)
				{
					tree[row + i + j].cost = tree[row + i].cost + tree[row + j].cost;
				}
			}
		}
		/*for (size_t s = 2; s < G; s++)
		{
			for (size_t s1 = s&(s-1); s1 > 1 ; s1 = s&(s1-1))
			{
				if (tree[row+s].cost > tree[row+s1].cost + tree[(s^s1)+row].cost )
				{
					tree[row+s].cost =  tree[row+s1].cost  + tree[(s^s1)+row].cost ;
					tree[row+s].p1 = s1;
					tree[row+s].p2 = s1^s;
				}
			}
		 }*/
	}
}
__global__ void one_label_lb(node *tree, size_t pitch_node, int *lb0, int N, int G, int inf, int *can_find_solution)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N)
	{
		node *row_u = (node *)((char *)tree + idx * pitch_node);
		for (size_t i = 1; i < G; i <<= 1)
		{
			if (lb0[idx] < row_u[i].cost)
			{
				lb0[idx] = row_u[i].cost;
			}
		}
		if (lb0[idx] != inf)
		{
			*can_find_solution = 1;
		}

		// for (size_t i = 1; i < G; i <<= 1)
		// {
		// 	for (size_t j = i << 1; j < G; j <<= 1)
		// 	{
		// 		if (row_u[i + j].cost > row_u[i].cost + row_u[j].cost)
		// 		{
		// 			row_u[i + j].cost = row_u[i].cost + row_u[j].cost;
		// 		}
		// 	}
		// }
		// for (size_t s = 2; s < G; s++)
		// {
		// 	for (size_t s1 = s&(s-1); s1 > 1 ; s1 = s&(s1-1))
		// 	{
		// 		if (row_u[s].cost > row_u[s1].cost + row_u[s^s1].cost)
		// 		{
		// 			row_u[s].cost = row_u[s1].cost + row_u[s^s1].cost;
		// 			row_u[s].p1 = s1;
		// 			row_u[s].p2 = s^s1;
		// 		}
		// 	}
		// }
	}
}
__global__ void dis0_init(queue_element *dis_queue, queue_element *new_dis_queue, node *tree, int *queue_size, int *new_queue_size, int *edge, int *edge_cost, int *pointer, size_t pitch_int, size_t pitch_node)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < *queue_size)
	{
		int u = dis_queue[idx].v, p = dis_queue[idx].p;
		node *row_u = (node *)((char *)tree + u * pitch_node);
		row_u[p].update = 0;
		for (int i = pointer[u]; i < pointer[u + 1]; i++)
		{
			int v = edge[i];
			int new_w = row_u[p].cost + edge_cost[i];
			node *row_v = (node *)((char *)tree + v * pitch_node);
			int old = atomicMin(&row_v[p].cost, new_w);
			if (new_w < old)
			{
				row_v[p].u = u;
					row_v[p].type = 1;
				int check = atomicCAS(&row_v[p].update, 0, 1);
				if (check == 0)
				{
					
					int pos = atomicAdd(new_queue_size, 1);
					new_dis_queue[pos] = {v, p};
				}
			}
		}
	}
}

__global__ void dis0_init_1(queue_element *dis_queue, queue_element *new_dis_queue, node *tree, int *queue_size, int *new_queue_size, int *edge, int *edge_cost, int *pointer, int width)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < *queue_size)
	{
		int u = dis_queue[idx].v, p = dis_queue[idx].p;
		tree[u * width + p].update = 0;
		for (int i = pointer[u]; i < pointer[u + 1]; i++)
		{
			int v = edge[i];
			int new_w = tree[u * width + p].cost + edge_cost[i];
			int old = atomicMin(&tree[v * width + p].cost, new_w);
			if (new_w < old)
			{
				tree[v * width + p].u = u;
				tree[v * width + p].type = 1;
				int check = atomicCAS(&tree[v * width + p].update, 0, 1);
				if (check == 0)
				{
					int pos = atomicAdd(new_queue_size, 1);
					new_dis_queue[pos] = {v, p};
				}
			}
		}
	}
}

graph_hash_of_mixed_weighted DPBF_GPU(node **host_tree, node *host_tree_one_d, CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph, int *pointer1, int *real_cost, int *community, int *c_size, non_overlapped_group_sets s, double *rt)
{
	// cudaSetDevice(1);
	int N = graph.V, width, height;
	int *queue_size, *new_queue_size, *best, *lb0, *non_overlapped_group_sets_IDs_gpu, *non_overlapped_group_sets_IDs_pointer_device, *can_find;
	int *all_pointer, *all_edge, *edge_cost;
	node *tree;
	queue_element *queue_device, *new_queue_device;
	queue_element *host_queue;
	double time_process = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	int G = cumpulsory_group_vertices.size();
	all_edge = graph.all_edge, all_pointer = graph.all_pointer, edge_cost = graph.all_edge_weight;
	int group_sets_ID_range = pow(2, G) - 1;
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	size_t avail(0); // 可用显存
	size_t total(0); // 总显存
	cudaMemGetInfo(&avail, &total);
	// std ::cout << "avail " << avail / 1024 / 1024 / 1024 << " total " << total / 1024 / 1024 / 1024 << std ::endl;
	//  printf("Number of devices: %d\n", deviceCount);

	long long unsigned int problem_size = N * pow(2, cumpulsory_group_vertices.size());

	cudaMallocManaged((void **)&queue_size, sizeof(int));
	cudaMallocManaged((void **)&new_queue_size, sizeof(int));
	cudaMallocManaged((void **)&can_find, sizeof(int));
	cudaMallocManaged((void **)&best, sizeof(int));
	cudaMallocManaged((void **)&lb0, N * sizeof(int));
	cudaMemset(lb0, 0, N * sizeof(int));
	cudaMallocManaged((void **)&queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&new_queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_gpu, sizeof(int) * s.length);
	cudaMalloc((void **)&non_overlapped_group_sets_IDs_pointer_device, sizeof(int) * (group_sets_ID_range + 3));
	cudaMemcpy(non_overlapped_group_sets_IDs_gpu, s.non_overlapped_group_sets_IDs.data(), sizeof(int) * s.length, cudaMemcpyHostToDevice);
	cudaMemcpy(non_overlapped_group_sets_IDs_pointer_device, s.non_overlapped_group_sets_IDs_pointer_host.data(), (group_sets_ID_range + 3) * sizeof(int), cudaMemcpyHostToDevice);

	width = group_sets_ID_range + 1, height = N;
	size_t pitch_node, pitch_int;
	host_queue = new queue_element[problem_size];
	cudaMallocPitch(&tree, &pitch_node, width * sizeof(node), height);
	// std::cout << "pitch " << pitch_node << " " << " width " << width << std::endl;
	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	time_process += runningtime;
	std ::cout << "main allocate cost time " << runningtime << std ::endl;
	*best = inf;
	begin = std::chrono::high_resolution_clock::now();
	/* 	for (size_t i = 0; i < G; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				host_dis[i][j] = inf;
			}
			// cout << "for group " << i << "  to ";
			for (size_t j = 0; j < group_graph[N + i].size(); j++)
			{

				cout << group_graph[N + i][j].first << " ";
				host_dis[i][group_graph[N + i][j].first] = 0;
			}
			// cout << "distance 0" << endl;
		}

		cudaMemcpy2D(dis, pitch_dis, host_dis, N * sizeof(int), N * sizeof(int), G, cudaMemcpyHostToDevice); */
	*queue_size = 0;
	set_max_ID(group_graph, cumpulsory_group_vertices, host_tree);
	for (int v = 0; v < N; v++)
	{
		host_tree[v][0].cost = 0;
		int group_set_ID_v = get_max(v, host_tree, width); /*time complexity: O(|Gamma|)*/
		for (size_t i = 1; i <= group_set_ID_v; i <<= 1)
		{
			if (i & group_set_ID_v)
			{

				host_tree[v][i].cost = 0;
				host_tree[v][i].type = 0;
				host_queue[*queue_size].v = v;
				host_queue[*queue_size].p = i;
				*queue_size += 1;
			}
		}
	}
	int *host_lb;
	host_lb = new int[N];
	cudaMemcpy2D(tree, pitch_node, host_tree_one_d, width * sizeof(node), width * sizeof(node), height, cudaMemcpyHostToDevice);
	*new_queue_size = 0;
	cudaMemcpy(queue_device, host_queue, *queue_size * sizeof(queue_element), cudaMemcpyHostToDevice);
	int r = 0, process = 0;

	while (*queue_size)
	{
		// std::cout << "main dis round " << r++ << " queue size " << *queue_size << std::endl;
		/* 		for (size_t i = 0; i < *queue_size; i++)
				{
					std::cout << " v " << queue_device[i].v << " p " << queue_device[i].p << "; ";
				}cout << endl; */

		dis0_init<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device, new_queue_device, tree, queue_size, new_queue_size, all_edge, edge_cost, all_pointer, pitch_int, pitch_node);
		cudaDeviceSynchronize();
		*queue_size = *new_queue_size;
		*new_queue_size = 0;
		std::swap(queue_device, new_queue_device);
	}
	*can_find = 0;

	one_label_lb<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, pitch_node, lb0, N, width, inf, can_find);
	cudaDeviceSynchronize();
	graph_hash_of_mixed_weighted solution_tree;
	if (*can_find == 0)
	{
		cout << "can not find a solution" << endl;
		return solution_tree;
	}
	cudaDeviceSynchronize();
	cudaMemcpy(host_lb, lb0, N * sizeof(int), cudaMemcpyDeviceToHost);

	// for (size_t i = 0; i < N; i++)
	// {
	// 	cout << i << " lb: " << lb0[i] << " ";
	// }
	// cudaMemcpy2D(host_dis, N * sizeof(int), dis, pitch_int, N * sizeof(int), G, cudaMemcpyDeviceToHost);
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	time_process += runningtime;
	std ::cout << "distance cost time " << runningtime << std ::endl;

	begin = std::chrono::high_resolution_clock::now();
	cudaMemcpy2D(host_tree_one_d, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
	*queue_size = 0, *new_queue_size = 0;

	for (int v = 0; v < N; v++)
	{

		host_tree[v][0].cost = 0;
		int group_set_ID_v = get_max(v, host_tree, width); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_sets_ID_range; p++)
		{ // p is non-empty; time complexity: O(2^|Gamma|) //get all its subset ,which is required in next merge and grow steps
			host_tree[v][p].update = 0;

			// host_tree[v][p].cost = inf;
			if ((p | group_set_ID_v) == group_set_ID_v)
			{ // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v
				/*T(v,p)*/
				host_tree[v][p].cost = 0;
				host_tree[v][p].type = 0;

				host_queue[*queue_size].v = v;
				host_queue[*queue_size].p = p;
				*queue_size += 1;
			}
		}
	}

	cudaMemcpy2D(tree, pitch_node, host_tree_one_d, width * sizeof(node), width * sizeof(node), height, cudaMemcpyHostToDevice);
	cudaMemcpy(queue_device, host_queue, *queue_size * sizeof(queue_element), cudaMemcpyHostToDevice);
	// std::cout << "queue size init " << *queue_size << std::endl;
	// std::cout << "queue init " << std::endl;
	/* 	for (size_t i = 0; i < *queue_size; i++)
		{
			std::cout << " v " << queue_device[i].v << " p " << queue_device[i].p << "; ";
		} */

	r = 0;
	begin = std::chrono::high_resolution_clock::now();

	while (*queue_size != 0)
	{
		process += *queue_size;
		std::cout << "round " << r++ << " queue size " << *queue_size << std::endl;
		//  cudaMemcpy(host_queue,queue_device,  *queue_size * sizeof(queue_element), cudaMemcpyDeviceToHost);
		//  		for (size_t i = 0; i < *queue_size; i++)
		//  		{
		//  			std::cout << " v " <<host_queue[i].v << " p " << host_queue[i].p << "; ";
		//  		}cout << endl;

		Relax<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device, queue_size, new_queue_device, new_queue_size, non_overlapped_group_sets_IDs_gpu,
																							 non_overlapped_group_sets_IDs_pointer_device, all_edge, edge_cost, all_pointer, pitch_node, tree, inf, best, group_sets_ID_range, lb0);
		// thrust::sort(queue_device, queue_device + *queue_size);
		cudaDeviceSynchronize();
		// cudaMemcpy2D(host_tree_one_d, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
		// 		for (size_t i = 0; i < N; i++)
		// 				{
		// 					cout << i << " ";
		// 					for (size_t j = 1; j <= group_sets_ID_range; j++)
		// 					{
		// 						cout << host_tree[i][j].cost << " ";
		// 					}
		// 					cout << endl;
		// 				}

		// cout << "new size = " << *new_queue_size << endl;
		*queue_size = *new_queue_size;
		*new_queue_size = 0;
		std::swap(queue_device, new_queue_device);
	}

	// std::cout << "while over,process node  " << process << std::endl;
	*pointer1 += process;
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	*rt = runningtime;
	std ::cout << "gpu cost time " << runningtime << std ::endl;
	cudaMemcpy2D(host_tree_one_d, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
	int min_cost = inf, min_node = -1;
	for (int i = 0; i < N; i++)
	{
		// cout << host_tree[i][group_sets_ID_range].cost << " ";
		if (host_tree[i][group_sets_ID_range].cost < min_cost)
		{
			min_cost = host_tree[i][group_sets_ID_range].cost;
			min_node = i;
		}
	}

	if (min_node == -1)
	{
		return solution_tree;
	}
	std::cout << "gpu root at " << min_node << " cost " << min_cost << std::endl;
	*real_cost = min_cost;
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
	cudaFree(queue_size);
	cudaFree(new_queue_size);
	cudaFree(non_overlapped_group_sets_IDs_gpu);
	cudaFree(non_overlapped_group_sets_IDs_pointer_device);
	cudaFree(tree);
	cudaFree(queue_device);
	cudaFree(new_queue_device);
	return solution_tree;
}
graph_hash_of_mixed_weighted DPBF_GPU_T(node **host_tree, node *host_tree_one_d, CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph, int *pointer1, int *real_cost, int *belong, std::vector<std::vector<int>> &community, non_overlapped_group_sets s, double *rt)
{
	// cudaSetDevice(1);
	auto pbegin = std::chrono::high_resolution_clock::now();
	auto begin = std::chrono::high_resolution_clock::now();
	auto pend = std::chrono::high_resolution_clock::now();
	int width, height, r = 0, process = 0, N = graph.V;
	int *queue_size, *new_queue_size, *best, *mid_queue_size, *lb0, *non_overlapped_group_sets_IDs_gpu, *non_overlapped_group_sets_IDs_pointer_device, *can_find;
	int *all_pointer, *all_edge, *edge_cost, *new_bfs_queue, *bfs_queue, *bfs_queue_size, *new_bfs_queue_size, *visited, mark_best;
	int *far_queue_size;
	node *tree;
	queue_element *queue_device, *new_queue_device, *mid_queue_device;
	queue_element *host_queue;
	double time_process = 0;
	int G = cumpulsory_group_vertices.size();
	all_edge = graph.all_edge, all_pointer = graph.all_pointer, edge_cost = graph.all_edge_weight;
	int group_sets_ID_range = pow(2, G) - 1;
	long long unsigned int problem_size = N * pow(2, cumpulsory_group_vertices.size());
	int *queue_size_p, *new_queue_size_p;
	queue_element *queue_device_p, *new_queue_device_p;
	cudaMallocManaged((void **)&queue_size_p, sizeof(int));
	cudaMallocManaged((void **)&new_queue_size_p, sizeof(int));
	cudaMallocManaged((void **)&far_queue_size, sizeof(int));
	cudaMallocManaged((void **)&can_find, sizeof(int));
	cudaMallocManaged((void **)&queue_device_p, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&new_queue_device_p, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&bfs_queue, N * sizeof(int));
	cudaMallocManaged((void **)&new_bfs_queue, N * sizeof(int));
	cudaMallocManaged((void **)&bfs_queue_size, sizeof(int));
	cudaMallocManaged((void **)&new_bfs_queue_size, sizeof(int));
	cudaMallocManaged((void **)&queue_size, sizeof(int));
	cudaMallocManaged((void **)&new_queue_size, sizeof(int));
	cudaMallocManaged((void **)&mid_queue_size, sizeof(int));
	cudaMallocManaged((void **)&best, sizeof(int));
	cudaMallocManaged((void **)&lb0, N * sizeof(int));
	cudaMallocManaged((void **)&visited, problem_size * sizeof(int));
	cudaMemset(lb0, 0, N * sizeof(int));
	cudaMemset(visited, 0, N * sizeof(int));
	cudaMallocManaged((void **)&queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&new_queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&mid_queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_gpu, sizeof(int) * s.length);
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_pointer_device, sizeof(int) * (group_sets_ID_range + 3));
	cudaMemcpy(non_overlapped_group_sets_IDs_gpu, s.non_overlapped_group_sets_IDs.data(), sizeof(int) * s.length, cudaMemcpyHostToDevice);
	cudaMemcpy(non_overlapped_group_sets_IDs_pointer_device, s.non_overlapped_group_sets_IDs_pointer_host.data(), (group_sets_ID_range + 3) * sizeof(int), cudaMemcpyHostToDevice);
	width = group_sets_ID_range + 1, height = N;
	host_queue = new queue_element[problem_size];
	cudaMallocManaged(&tree, width * height * sizeof(node));

	// std::cout << "pitch " << pitch_node << " " << " width " << width << std::endl;
	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	time_process += runningtime;
	std ::cout << "main allocate cost time " << runningtime << std ::endl;
	size_t avail(0); // 可用显存
	size_t total(0); // 总显存
	cudaMemGetInfo(&avail, &total);
	std ::cout << "avail " << avail / 1024 / 1024 / 1024 << " total " << total / 1024 / 1024 / 1024 << std ::endl;
	*best = inf;
	begin = std::chrono::high_resolution_clock::now();
	*queue_size = 0;

	/*
	int minimum_group = *cumpulsory_group_vertices.begin();
	for (const auto &key : cumpulsory_group_vertices)
	{
		if (group_graph[key].size() < group_graph[minimum_group].size())
		{
			minimum_group = key;
		}
	}
	*bfs_queue_size = 0;
	*queue_size_p = 0;
	*new_bfs_queue_size = 0;
	for (size_t i = 0; i < group_graph[minimum_group].size(); i++)
	{
		bfs_queue[*bfs_queue_size] = group_graph[minimum_group][i].first;
		*bfs_queue_size += 1;
	}

	for (size_t i = 0; i < 1; i++)
	{
		std ::cout << std ::endl;
		BFS<<<(*bfs_queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(bfs_queue, bfs_queue_size, new_bfs_queue, new_bfs_queue_size, all_edge, visited, all_pointer);
		cudaDeviceSynchronize();
		*bfs_queue_size = *new_bfs_queue_size;
		*new_bfs_queue_size = 0;
		std::swap(bfs_queue, new_bfs_queue);
	}*/
	set_max_ID(group_graph, cumpulsory_group_vertices, host_tree);

	for (int v = 0; v < N; v++)
	{
		host_tree[v][0].cost = 0;
		int group_set_ID_v = get_max(v, host_tree, width); /*time complexity: O(|Gamma|)*/
		for (size_t i = 1; i <= group_set_ID_v; i <<= 1)
		{
			if (i & group_set_ID_v)
			{
				host_tree[v][i].cost = 0;
				host_tree[v][i].type = 0;
				host_queue[*queue_size].v = v;
				host_queue[*queue_size].p = i;
				*queue_size += 1;
			}
		}
	}
	int *host_lb;
	host_lb = new int[N];
	cudaMemcpy(tree, host_tree_one_d, width * sizeof(node) * height, cudaMemcpyHostToDevice);
	*new_queue_size = 0;
	cudaMemcpy(queue_device, host_queue, *queue_size * sizeof(queue_element), cudaMemcpyHostToDevice);

	while (*queue_size)
	{
		// std::cout << "dis round " << r++ << " queue size " << *queue_size << std::endl;
		dis0_init_1<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device, new_queue_device, tree, queue_size, new_queue_size, all_edge, edge_cost, all_pointer, width);
		cudaDeviceSynchronize();
		*queue_size = *new_queue_size;
		*new_queue_size = 0;
		std::swap(queue_device, new_queue_device);
	}
	*can_find = 0;
	one_label_lb_1<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(tree, width, lb0, N, width, inf, can_find);
	cudaDeviceSynchronize();
	graph_hash_of_mixed_weighted solution_tree;
	if (*can_find == 0)
	{
		cout << "can not find a solution" << endl;
		return solution_tree;
	}
	cudaMemcpy(host_lb, lb0, N * sizeof(int), cudaMemcpyDeviceToHost);
	// for (size_t i = 0; i < N; i++)
	// {
	// 	cout << i << " lb: " << lb0[i] << " ";
	// }
	// cout<<endl;
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	time_process += runningtime;
	std ::cout << "distance cost time " << runningtime << std ::endl;
	begin = std::chrono::high_resolution_clock::now();
	cudaMemcpy(host_tree_one_d, tree, width * sizeof(node) * height, cudaMemcpyDeviceToHost);
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	std ::cout << "till copy1 cost time " << runningtime << std ::endl;
	*queue_size = 0, *new_queue_size = 0;
	for (int v = 0; v < N; v++)
	{
		// int dev = belong[v];
		host_tree[v][0].cost = 0;

		int group_set_ID_v = get_max(v, host_tree, width); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_sets_ID_range; p++)
		{ // p is non-empty; time complexity: O(2^|Gamma|) //get all its subset ,which is required in next merge and grow steps
			host_tree[v][p].update = 0;

			// host_tree[v][p].cost = inf;
			if ((p | group_set_ID_v) == group_set_ID_v)
			{ // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v
				/*T(v,p)*/
				host_tree[v][p].cost = 0;
				host_tree[v][p].type = 0;

				// queue_device_p[queue_size_p[dev]].v = v;
				// queue_device_p[queue_size_p[dev]].p = p;
				// queue_size_p[dev] += 1;

				if (visited[v] == 1)
				{
					queue_device_p[*queue_size_p].v = v;
					queue_device_p[*queue_size_p].p = p;
					*queue_size_p += 1;
					queue_device[*queue_size].v = v;
					queue_device[*queue_size].p = p;
					*queue_size += 1;
				}
				else
				{
					queue_device[*queue_size].v = v;
					queue_device[*queue_size].p = p;
					*queue_size += 1;
				}
			}
		}
	}
end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	std ::cout << "till mark cost time " << runningtime << std ::endl;
	cudaMemcpy(tree, host_tree_one_d, width * sizeof(node) * height, cudaMemcpyHostToDevice);
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	std ::cout << "till copy2cost time " << runningtime << std ::endl;
	// std::cout << "queue size init " << *queue_size << std::endl;
	// std::cout << "queue init " << std::endl;
	/* 	for (size_t i = 0; i < *queue_size; i++)
		{
			std::cout << " v " << queue_device[i].v << " p " << queue_device[i].p << "; ";
		} */

	r = 0;
	int up = 125932;
	double percent = 0.25;

	
	// || *queue_size_p != 0
	while (*queue_size != 0 && r < 20)
	{
		process += *queue_size;
		 std::cout << "round " << r++ << " queue size " << *queue_size << std::endl;
		/*cudaMemcpy(host_queue,queue_device,  *queue_size * sizeof(queue_element), cudaMemcpyDeviceToHost);
				for (size_t i = 0; i < *queue_size; i++)
				{
					std::cout << " v " <<host_queue[i].v << " p " << host_queue[i].p << "; ";
				}cout << endl;*/

		Relax_1<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device, queue_size, new_queue_device, new_queue_size, non_overlapped_group_sets_IDs_gpu,
																							   non_overlapped_group_sets_IDs_pointer_device, all_edge, edge_cost, all_pointer, width, tree, inf, best, group_sets_ID_range, lb0);

		// if (*queue_size_p)
		// {
		// 	Relax_1<<<(*queue_size_p + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device_p, queue_size_p, new_queue_device_p, new_queue_size_p, non_overlapped_group_sets_IDs_gpu,
		// 																							 non_overlapped_group_sets_IDs_pointer_device, all_edge, edge_cost, all_pointer, width, tree, inf, best, group_sets_ID_range, lb0);
		// }

		cudaDeviceSynchronize();
		cudaMemset(visited, 0, *new_queue_size * sizeof(int));
		if (*new_queue_size > up)
		{ // 如果新队列的长度大于运行能力 则运行近似优先
			//	cout << *new_queue_size << "> up" << endl;
			*queue_size = 0;
			pre_scan<<<(*new_queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(new_queue_device, new_queue_size, queue_device, queue_size, best, up, percent, visited, mid_queue_device, mid_queue_size);
			cudaDeviceSynchronize();
			if (*queue_size < up)
			{
				cout << *queue_size << "< up need fill up , percent = " << percent << "best = " << *best << endl;
				// 经过pre_scan得到的queuesize<up 说明筛出来的太少了 要从远堆中补充得到满队 下次要增大门槛
				int need = up - *queue_size;
				cudaMemcpy(queue_device + *queue_size, mid_queue_device + *mid_queue_size - need, need * sizeof(queue_element), cudaMemcpyDeviceToDevice);
				*mid_queue_size -= need;
				percent = 0.5*double(need)/up;
				*queue_size = up;
				
				// percent=min(percent+double(need)/double(up),1.0);
			}
			else
			{ // 经过pre_scan得到的queuesize>up 说明筛出来的太多了 下次要减少门槛
				cout << *queue_size << "> up need cut dow, percent = " << percent << " best = " << *best << endl;
				double cut = *queue_size-up;
				percent = 0.5*cut/(*queue_size);
				*queue_size = up;
				// cout << *queue_size << ">up" << endl;
				
			}

			/*	cudaMemcpy2D(host_tree_one_d, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
						for (size_t i = 0; i < N; i++)
								{
									cout << i << " ";
									for (size_t j = 1; j <= group_sets_ID_range; j++)
									{
										cout << host_tree[i][j].cost << " ";
									}
									cout << endl;
								}*/

			// cout << "new size = " << *new_queue_size << endl;

			*new_queue_size = *mid_queue_size;
			std::swap(new_queue_device, mid_queue_device);
			*mid_queue_size = 0;
			// *queue_size_p = *new_queue_size_p;
			// *new_queue_size_p = 0;
			// std::swap(queue_device_p, new_queue_device_p);
			// if (*queue_size_p == 0)
			// {
			// 	pend = std::chrono::high_resolution_clock::now();
			// 	mark_best = *best;
			// 	*rt = runningtime;
			// }
			cout << "near size = " << *queue_size << endl;
			cout << "far size = " << *new_queue_size << endl;
		}
		else
		{
			*queue_size = *new_queue_size;
			*new_queue_size = 0;
			std::swap(queue_device, new_queue_device);
		}
	}

	// std::cout << "while over,process node  " << process << std::endl;
	*pointer1 += process;
	end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	std ::cout << "gpu cost time " << runningtime << std ::endl;
	begin = std::chrono::high_resolution_clock::now();
	cudaMemcpy(host_tree_one_d, tree, width * sizeof(node) * height, cudaMemcpyDeviceToHost);
	int min_cost = inf, min_node = -1;
	for (int i = 0; i < N; i++)
	{
		// cout << host_tree[i][group_sets_ID_range].cost << " ";
		if (host_tree[i][group_sets_ID_range].cost < min_cost)
		{
			min_cost = host_tree[i][group_sets_ID_range].cost;
			min_node = i;
		}
	}
	*real_cost = min_cost;
	std::cout << "gpu root at " << min_node << " cost " << min_cost << std::endl;
	if (min_node == -1)
	{
		return solution_tree;
	}

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
	if (min_cost == mark_best)
	{
		runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(pend - pbegin).count() / 1e9; // s
		cout << "part find solution" << endl;
		*rt = runningtime;
	}
	else
	{
		runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - pbegin).count() / 1e9; // s
		*rt = runningtime;
	}
	cudaFree(tree);
	cudaFree(queue_device);
	cudaFree(new_queue_device);
	cudaFree(mid_queue_device);
	cudaFree(bfs_queue);
	cudaFree(new_bfs_queue);
		end = std::chrono::high_resolution_clock::now();
	runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	std ::cout << "form tree cost time " << runningtime << std ::endl;
	return solution_tree;
}
