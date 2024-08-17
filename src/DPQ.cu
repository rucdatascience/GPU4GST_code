
#include <DPQ.cuh>
#include <thrust/device_vector.h>
using namespace std;
/*this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.

time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )*/
typedef struct queue_element
{
	int v, p;
} queue_element;
typedef struct node
{
	int update = 0;
	int type;	  // =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge
	int cost, lb; // cost of this tree T(v,p);
	int u;		  // if this tree is built by grown, then it's built by growing edge (v,u);
	int p1, p2;	  // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);
} node;
int E, N, width, height;
int *lb1, *lb2;
int *visit, *queue_size, *tree_cost, *new_queue_size, *best, *dis, *in_queue_check, *dis_queue, *new_dis_queue;
int *all_pointer, *all_edge, *edge_cost, *non_overlapped_group_sets_IDs_gpu, *non_overlapped_group_sets_IDs_pointer_device;
dim3 blockPerGrid, threadPerGrid;
node *tree;
queue_element *queue_device, *new_queue_device;
std::vector<int> non_overlapped_group_sets_IDs_pointer_host;
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

void graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(int group_sets_ID_range)
{

	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID;

	time complexity: O(4^|Gamma|), since group_sets_ID_range=2^|Gamma|;

	the original DPBF code use the same method in this function, and thus has the same O(4^|Gamma|) complexity;*/

	std::vector<int> non_overlapped_group_sets_IDs; // <set_ID, non_overlapped_group_sets_IDs>
	int len = 0;
	for (int i = 1; i <= group_sets_ID_range; i++)
	{ // i is a nonempty group_set ID
		non_overlapped_group_sets_IDs_pointer_host[i] = len;
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
	non_overlapped_group_sets_IDs_pointer_host[group_sets_ID_range + 1] = len;

	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_gpu, sizeof(int) * len);
	cudaMemcpy(non_overlapped_group_sets_IDs_gpu, non_overlapped_group_sets_IDs.data(), sizeof(int) * len, cudaMemcpyHostToDevice);
	std::cout << "len= " << len << std::endl;
}
__device__ void get_lb(int *edge, int *edge_cost, int *pointer, size_t pitch_node, size_t pitch_int, size_t pitch_dis, int *dis, node *tree, int *tree_cost, int inf, int *best, int full, int cost, int *lb, int v, int p)
{
	int *row_v = (int *)((char *)tree_cost + v * pitch_int);
	node *row_node_v = (node *)((char *)tree + v * pitch_node);

	int x_slash = full - p, g = 0, one_label = -1;
	while (x_slash)
	{
		if (x_slash & 1)
		{
			int *row_g = (int *)((char *)dis + g * pitch_int);
			if (row_g[v] > one_label)
			{
				one_label = row_g[v];
			}
		}

		g++;
		x_slash >>= 1;
	}
	*lb = one_label + cost;
}

__global__ void Relax(queue_element *Queue_dev, int *queue_size, queue_element *new_queue_device, int *new_queue_size, int *sets_IDs, int *sets_IDS_pointer, int *edge, int *edge_cost, int *pointer, size_t pitch_node, size_t pitch_int, size_t pitch_dis, int *dis, node *tree, int *tree_cost, int inf, int *best, int full,int *lb1,int *lb2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < *queue_size)
	{
		queue_element top_node = Queue_dev[idx];
		int v = top_node.v, p = top_node.p;
		int *row_v = (int *)((char *)tree_cost + v * pitch_int);
		node *row_node_v = (node *)((char *)tree + v * pitch_node);
		int x_slash = full - p;
		if (row_v[x_slash] != inf && row_v[p] != inf)
		{
			int new_best = row_v[x_slash] + row_v[p];
			atomicMin(best, new_best);
			atomicMin(&row_v[full], new_best);
			atomicMin(&row_node_v[full].cost, new_best);
			int check = atomicCAS(&row_node_v[full].update, 0, 1);
			if (!check)
			{
				int pos = atomicAdd(new_queue_size, 1);
				new_queue_device[pos] = {v, full};
			}
		}

		if (row_v[p] > (*best) / 2)
		{
			return;
		}
		for (int i = pointer[v]; i < pointer[v + 1]; i++)
		{
			/*grow*/
			int lb;
			int u = edge[i];
			int cost_euv = edge_cost[i];
			int *row_u = (int *)((char *)tree_cost + u * pitch_int);
			int *row_lb1 = (int*)((char *)lb1 + u * pitch_int);
			int *row_lb2 = (int*)((char *)lb2 + u * pitch_int);
			int grow_tree_cost = row_v[p] + cost_euv;
			int old = atomicMin(&row_u[p], grow_tree_cost);
			node *row_node_u = (node *)((char *)tree + u * pitch_node);
			// get_lb(edge,edge_cost,pointer,pitch_node,pitch_int,pitch_dis,dis,tree,tree_cost,inf,best,full,grow_tree_cost,plb,v,p);
			x_slash = full - p;
			lb = row_lb1[x_slash]>row_lb2[x_slash]?row_lb1[x_slash]:row_lb2[x_slash];
			int g = 0, one_label = -1;
			if (p == full)
			{
				atomicMin(best, grow_tree_cost);
			}
			while (x_slash)
			{
				if (x_slash & 1)
				{
					int *row_g = (int *)((char *)dis + g * pitch_dis);
					if (row_g[u] > one_label)
					{
						one_label = row_g[u];
					}
				}
				g++;
				x_slash = x_slash >> 1;
			}
			lb = lb>one_label?lb:one_label;

			atomicMin(&row_node_u[p].cost, grow_tree_cost);
			int low_bound = atomicMax(&row_node_u[p].lb, lb);

			if (low_bound > *best)
			{
				continue;
			}

			if (old >= grow_tree_cost && grow_tree_cost != inf)
			{
				row_node_u[p].type = 1;
				row_node_u[p].u = v;
				row_node_u[p].cost = grow_tree_cost;
				// enqueue operation
				int check = atomicCAS(&row_node_u[p].update, 0, 1);
				if (!check)
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
			int *row_v = (int *)((char *)tree_cost + v * pitch_int);
			int cost_Tvp1 = row_v[p1], cost_Tvp2 = row_v[p2];
			int p1_cup_p2 = p1 + p2;
			int merged_tree_cost = cost_Tvp1 + cost_Tvp2;
			// && merged_tree_cost < 2 / 3 * (*best)
			int old = atomicMin(&row_v[p1_cup_p2], merged_tree_cost);
			atomicMin(&row_node_v[p1_cup_p2].cost, merged_tree_cost);
			int *row_lb1 = (int*)((char *)lb1 + v * pitch_int);
			int *row_lb2 = (int*)((char *)lb2 + v * pitch_int);
			x_slash = full - p;
			int lb;
			lb = row_lb1[x_slash]>row_lb2[x_slash]?row_lb1[x_slash]:row_lb2[x_slash];
			int g = 0, one_label = -1;

			while (x_slash)
			{
				if (x_slash & 1)
				{
					int *row_g = (int *)((char *)dis + g * pitch_dis);
					if (row_g[v] > one_label)
					{
						one_label = row_g[v];
					}
				}
				g++;
				x_slash = x_slash >> 1;
			}
			lb = lb>one_label?lb:one_label;
			int low_bound = atomicMax(&row_node_v[p].lb, lb);

			if (low_bound > *best)
			{
				continue;
			}
			if (old >= merged_tree_cost && merged_tree_cost != inf)
			{ // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				row_node_v[p1_cup_p2].type = 2;
				row_node_v[p1_cup_p2].p1 = p1;
				row_node_v[p1_cup_p2].p2 = p2;
				row_node_v[p1_cup_p2].cost = merged_tree_cost;

				if (merged_tree_cost < 0.667 * (*best))
				{
					int check = atomicCAS(&row_node_v[p1_cup_p2].update, 0, 1);
					if (!check)
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

__global__ void dis_init(int v, int *in_queue_check, int pitch_int, int N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N)
	{
		int *row_v = (int *)((char *)in_queue_check + v * pitch_int);
		row_v[idx] = 0;
	}
}
__global__ void dis_Relax(int *dis_queue, int *new_dis_queue, int *dis, int *in_queue_check, int *queue_size, int *new_queue_size, int *edge, int *edge_cost, int *pointer, size_t pitch_int, int source)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < *queue_size)
	{
		int u = dis_queue[idx];
		int *row_s = (int *)((char *)dis + source * pitch_int);
		int *row_in = (int *)((char *)in_queue_check + source * pitch_int);

		for (int i = pointer[u]; i < pointer[u + 1]; i++)
		{
			int v = edge[i];
			int new_w = row_s[u] + edge_cost[i];

			int old = atomicMin(&row_s[v], new_w);
			if (new_w < old)
			{
				int check = atomicCAS(&row_in[v], 0, 1);
				if (!check)
				{
					int pos = atomicAdd(new_queue_size, 1);
					new_dis_queue[pos] = v;
				}
			}
		}
	}
}
graph_hash_of_mixed_weighted DPBF_GPU(CSR_graph &graph, std::unordered_set<int> &cumpulsory_group_vertices, graph_v_of_v_idealID &group_graph, graph_v_of_v_idealID &input_graph)
{

	E = graph.E_all, N = graph.V;
	int G = cumpulsory_group_vertices.size();
	all_edge = graph.all_edge, all_pointer = graph.all_pointer, edge_cost = graph.all_edge_weight;
	int group_sets_ID_range = pow(2, G) - 1;
	int inf = 1024;
	non_overlapped_group_sets_IDs_pointer_host.resize(group_sets_ID_range + 3);
	long long unsigned int problem_size = N * pow(2, cumpulsory_group_vertices.size());
	cudaMallocManaged((void **)&non_overlapped_group_sets_IDs_pointer_device, sizeof(int) * (group_sets_ID_range + 3));
	cudaMallocManaged((void **)&queue_size, sizeof(int));
	cudaMallocManaged((void **)&new_queue_size, sizeof(int));
	cudaMallocManaged((void **)&best, sizeof(int));
	cudaMallocManaged((void **)&queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&new_queue_device, problem_size * sizeof(queue_element));
	cudaMallocManaged((void **)&dis_queue, problem_size * sizeof(int));
	cudaMallocManaged((void **)&in_queue_check, problem_size * sizeof(int));
	cudaMallocManaged((void **)&new_dis_queue, problem_size * sizeof(int));
	graph_v_of_v_idealID_DPBF_non_overlapped_group_sets_gpu(group_sets_ID_range);
	cudaMemcpy(non_overlapped_group_sets_IDs_pointer_device, non_overlapped_group_sets_IDs_pointer_host.data(), (group_sets_ID_range + 3) * sizeof(int), cudaMemcpyHostToDevice);
	std::cout << "group range " << group_sets_ID_range << std::endl;
	threadPerGrid.x = THREAD_PER_BLOCK;
	blockPerGrid.x = (N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK;
	width = group_sets_ID_range + 1, height = N;
	size_t pitch_node, pitch_int, pitch_vis, pitch_dis;
	node host_tree[height][width];
	int host_cost[height][width];
	int host_dis[G][N], host_f1[N][width], host_f2[N][width];
	queue_element host_queue[problem_size];
	cudaMallocPitch(&dis, &pitch_dis, N * sizeof(int), G);
	cudaMallocPitch(&in_queue_check, &pitch_vis, N * sizeof(int), G);
	cudaMallocPitch(&tree, &pitch_node, width * sizeof(node), height);
	cudaMallocPitch(&tree_cost, &pitch_int, width * sizeof(int), height);
	cudaMallocPitch(&lb1, &pitch_int, width * sizeof(int), G);
	cudaMallocPitch(&lb2, &pitch_int, width * sizeof(int), G);
	// cudaMemset3D(devPitchedPtr, inf, extent);
	std::cout << "pitch " << pitch_node << " " << " width " << width << std::endl;

	*best = inf;
	for (size_t i = 0; i < G; i++)
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

	cudaMemcpy2D(dis, pitch_dis, host_dis, N * sizeof(int), N * sizeof(int), G, cudaMemcpyHostToDevice);
	for (size_t i = 0; i < G; i++)
	{

		// cout << "queue init for " << i << " group" << endl;
		*queue_size = group_graph[N + i].size();
		*new_queue_size = 0;
		for (size_t j = 0; j < group_graph[N + i].size(); j++)
		{
			dis_queue[j] = group_graph[N + i][j].first;
			// cout << dis_queue[j] << " ";
		}
		dis_init<<<(N + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(i, in_queue_check, pitch_vis, N);
		while (*queue_size)
		{

			// cout << "size " << *queue_size << endl;
			dis_Relax<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(dis_queue, new_dis_queue, dis, in_queue_check, queue_size, new_queue_size, all_edge, edge_cost, all_pointer, pitch_dis, i);
			cudaDeviceSynchronize();
			/* 			cudaMemcpy2D(host_dis, N * sizeof(int), dis, pitch_dis, N * sizeof(int), G, cudaMemcpyDeviceToHost);
						cudaDeviceSynchronize();
						for (size_t j = 0; j < N; j++)
						{
							cout << j << ":" << host_dis[i][j] << " ";
						}
						cout << endl; */

			*queue_size = *new_queue_size;
			*new_queue_size = 0;
			std::swap(dis_queue, new_dis_queue);
		}
	}
	cudaMemcpy2D(host_dis, N * sizeof(int), dis, pitch_int, N * sizeof(int), G, cudaMemcpyDeviceToHost);
	int host_w[G][G][width], host_w1[G][width];
	for (size_t i = 0; i < G; i++)
	{
		for (size_t j = 0; j < G; j++)
		{
			for (size_t k = 0; k <= width; k++)
			{

				host_w[i][j][k] = inf;
			}
		}
	}
	for (size_t vi = 0; vi < G; vi++)
	{
		for (size_t j = 0; j < group_graph[vi].size(); j++)
		{
			for (size_t vg = 0; vg < G; vg++)
			{
				for (size_t ii = 0; ii < group_graph[vg].size(); ii++)
				{
					host_w[vi][vg][0] = min(host_w[vi][vg][0], host_dis[group_graph[vi][j].first][group_graph[vg][ii].first]);
				}
			}
		}
	}
	for (size_t i = 0; i < G; i++)
	{
		for (size_t j = 0; j < G; j++)
		{
			if (i == j)
			{
				for (size_t x = 1; x <= group_sets_ID_range; x++)
				{
					host_w[i][j][x] = 0;
				}
			}
			else
			{
				for (size_t x = 1; x <= group_sets_ID_range; x++)
				{
					int p_except = 1;
					while (p_except <= x)
					{
						if (p_except & x == p_except)
						{
							host_w[i][j][x] = min(host_w[i][j][x], host_dis[i][p_except] + host_w[i][p_except][p_except ^ x]);
						}
						p_except <<= 1;
					}
				}
			}
		}
	}
	for (size_t v = 0; v < N; v++)
	{
		for (size_t x = 1; x < width; x++)
		{
			int vi = 1;
			while (vi <= x)
			{
				int vj = vi << 1;
				while (vj <= x)
				{
					host_f1[v][x] = min(host_f1[v][x], host_dis[vi][v] + host_dis[vj][v] + host_w[vi][vj][x]);
					vj <<= 1;
				}

				vi <<= 1;
			}
			host_f1[v][x] /= 2;
		}
	}
	for (size_t i = 0; i < G; i++)
	{
		for (size_t x = 1; x < width; x++)
		{
			for (size_t k = 0; k < G; k++)
			{
				host_w1[i][x] = min(host_w1[i][x], host_w[i][k][x]);
			}
		}
	}

	for (size_t v = 0; v < N; v++)
	{
		for (size_t x = 1; x < width; x++)
		{
			int vj = 1;
			int dist = inf;
			while (vj <= x)
			{
				if (vj & x == vj)
				{
					dist = min(dist, host_dis[vj][v]);
				}
				vj <<= 1;
			}
			int vi = 1;
			while (vi <= x)
			{
				if (vi & x == vi)
				{
					host_f2[v][x] = max(host_f2[v][x], host_dis[vi][v] + host_w1[vi][x] + dist);
				}
				vi <<= 1;
			}
			host_f2[v][x] /= 2;
		}
	}
	cudaMemcpy2D(lb1, pitch_int, host_f1, width * sizeof(int), width * sizeof(int), N, cudaMemcpyHostToDevice);
	cudaMemcpy2D(lb2, pitch_int, host_f2, width * sizeof(int), width * sizeof(int), N, cudaMemcpyHostToDevice);
	/* 	cudaExtent extent = make_cudaExtent(N * sizeof(int), N, width);
		cudaPitchedPtr devPitchedPtr;
		cudaMalloc3D(&devPitchedPtr, extent);
		cudaMemcpy3DParms DevToHost = {0};
		DevToHost.srcPtr = devPitchedPtr;
		DevToHost.dstPtr = make_cudaPitchedPtr((void *)host_w, width * sizeof(int), width, N);
		DevToHost.extent = extent;
		DevToHost.kind = cudaMemcpyDeviceToHost;
		cudaMemcpy3D(&DevToHost); */

	/* 	for (size_t i = 0; i < G; i++)
		{
			cout << i << " ";
			for (size_t j = 0; j < N; j++)
			{
				cout << j << ":" << host_dis[i][j] << " ";
			}
			cout << endl;
		} */

	*queue_size = 0, *new_queue_size = 0;
	for (int v = 0; v < N; v++)
	{
		host_tree[v][0].cost = 0;
		host_cost[v][0] = 0;
		int group_set_ID_v = graph_v_of_v_idealID_DPBF_vertex_group_set_ID_gpu(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_sets_ID_range; p++)
		{ // p is non-empty; time complexity: O(2^|Gamma|) //get all its subset ,which is required in next merge and grow steps
			host_tree[v][p].cost = inf;
			host_cost[v][p] = inf;
			if ((p | group_set_ID_v) == group_set_ID_v)
			{ // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v
				/*T(v,p)*/
				host_tree[v][p].cost = 0;
				host_tree[v][p].type = 0;
				host_cost[v][p] = 0;
				host_queue[*queue_size].v = v;
				host_queue[*queue_size].p = p;
				*queue_size += 1;
			}
		}
	}

	cudaMemcpy2D(tree, pitch_node, host_tree, width * sizeof(node), width * sizeof(node), height, cudaMemcpyHostToDevice);
	cudaMemcpy2D(tree_cost, pitch_int, host_cost, width * sizeof(int), width * sizeof(int), height, cudaMemcpyHostToDevice);
	cudaMemcpy(queue_device, host_queue, *queue_size * sizeof(queue_element), cudaMemcpyHostToDevice);
	std::cout << "queue size init " << *queue_size << std::endl;
	std::cout << "queue init " << std::endl;
	for (size_t i = 0; i < *queue_size; i++)
	{
		std::cout << " v " << queue_device[i].v << " p " << queue_device[i].p << "; ";
	}
	cout << endl;
	int r = 0;
	while (*queue_size != 0)
	{
		std::cout << "round " << r++ << std::endl;
		std::cout << "queue size " << *queue_size << std::endl;
		/* 		for (size_t i = 0; i < *queue_size; i++)
				{
					std::cout << " v " << queue_device[i].v << " p " << queue_device[i].p << "; ";
				} */
		cout << endl;
		Relax<<<(*queue_size + THREAD_PER_BLOCK - 1) / THREAD_PER_BLOCK, THREAD_PER_BLOCK>>>(queue_device, queue_size, new_queue_device, new_queue_size, non_overlapped_group_sets_IDs_gpu,
																							 non_overlapped_group_sets_IDs_pointer_device, all_edge, edge_cost, all_pointer, pitch_node, pitch_int, pitch_dis, dis, tree, tree_cost, inf, best, group_sets_ID_range,lb1,lb2);
		cudaDeviceSynchronize();
		/* 		cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
				cudaMemcpy2D(host_cost, width * sizeof(int), tree_cost, pitch_int, width * sizeof(int), height, cudaMemcpyDeviceToHost);
				for (size_t i = 0; i < N; i++)
				{
					cout << i << " ";
					for (size_t j = 1; j <= group_sets_ID_range; j++)
					{
						cout << host_cost[i][j] << " ";
					}
					cout << endl;
				}
				cout<<"new size = "<<*new_queue_size<<endl; */
		*queue_size = *new_queue_size;
		*new_queue_size = 0;
		std::swap(queue_device, new_queue_device);
	}
	std::cout << "while over" << std::endl;
	cudaMemcpy2D(host_tree, width * sizeof(node), tree, pitch_node, width * sizeof(node), height, cudaMemcpyDeviceToHost);
	cudaMemcpy2D(host_cost, width * sizeof(int), tree_cost, pitch_int, width * sizeof(int), height, cudaMemcpyDeviceToHost);

	std::cout << "all copy complete ,now list cost " << std::endl;
	int min_cost = inf, min_node = -1;
	for (int i = 0; i < N; i++)
	{
		// cout << host_tree[i][group_sets_ID_range].cost << " ";
		if (host_cost[i][group_sets_ID_range] < min_cost)
		{
			min_cost = host_cost[i][group_sets_ID_range];
			min_node = i;
		}
	}

	std::cout << "root " << min_node << "cost " << min_cost << std::endl;
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
