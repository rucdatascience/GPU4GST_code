# rucgraph

This is a library of cpp header files mainly used for performing graph computing tasks. Here is a [Table of Contents](assets/Introduction.pdf). Some other descriptions can be found on this [page](https://yahuisun.github.io/rucgraph/). <b>Welcome to joining the development of rucgraph by pulling requests!</b> (Contact: yahuisun@ruc.edu.cn)



<p align="center">
<img src="/assets/images/manual_Kaggle_lowquality.jpg" alt="drawing" width="200"/>
</p>



## How to use rucgraph

<b>Normally, there is an example in each header file to explain how to use the codes in this file. If there is not, you can add a mark, e.g., "/\*an example is needed\*/", in the end of this file, and then pull a request. I will add an example ASAP.</b>

You can download the whole repository into a rucgraph folder, and add the path of this folder when compiling cpp codes.


### An example of using rucgraph on a Linux server

1. Download the whole rucgraph folder, and save it on the server, e.g., the path of the saved folder is "root/rucgraph".

2. Given a cpp file named as "try.cpp", the contents in which are
```
using namespace std;

#include <data_structures/PairingHeapYS_with_offset.h>

int main()
{
	example_PairingHeapYS_with_offset();
}
```
, where "PairingHeapYS_with_offset.h" is a header file in rucgraph. In the terminal, compile and run the above codes using the following commands:
```
g++ -std=c++17 -I/root/rucgraph try.cpp a.out
./a.out
```
, where "-I/root/rucgraph" is to add the path of the rucgraph folder when compiling. Then, you successfully run an example application of an augmented pairing heap, detailed in "PairingHeapYS_with_offset.h".
![](/assets/images/202212171254231.png)


### An example of using both boost and rucgraph on a Linux server

Some header files in rucgraph require [the Boost library](https://www.boost.org/). To use these files, it is required to download the Boost library as well. For example, given a cpp file named as "try.cpp", the contents in which are
```
#include <graph_v_of_v_idealID/common_algorithms/graph_v_of_v_idealID_shortest_paths.h>

int main()
{
	test_graph_v_of_v_idealID_shortest_paths();
}
```
In the terminal, compile and run the above codes using the following commands:
```
g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -o a.out
./a.out
```
, where "-I/home/boost_1_75_0" is to add the path of the boost folder when compiling. Then, you successfully run a test of finding shortest paths in a graph, detailed in "graph_v_of_v_idealID_shortest_paths.h".
![](/assets/images/20221217183837.png)

<br/>


# Highlights

## A Pairing Heap that can change all keys in O(1) time

rucgraph contains <span style="color:red">an augmented Pairing Heap that can change all keys in O(1) time!</span> It achieves this by associating an offset value with each node in the heap. To use this heap, include the following header file in your codes.
```
#include <data_structures/PairingHeapYS_with_offset.h>
```
<br/>


# Notes for updating files.

1. It is preferable to add an example in the front or end of each h file, for describing how to use the codes in the file.
2. The recorded codes, except those in the build_in_progress folder, should be as stable as possible.
