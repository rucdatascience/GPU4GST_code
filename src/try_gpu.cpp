#include <iostream>
#include <fstream>

using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
#include "test_gpu.h"
int main()
{
	test_graph_v_of_v_idealID_DPBF_only_ec_gpu();
}