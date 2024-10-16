#pragma once

#include <text_mining/parse_string.h>

void graph_hash_of_mixed_weighted_read_for_GSTP(std::string instance_name,
												graph_hash_of_mixed_weighted &input_graph, graph_hash_of_mixed_weighted &group_graph,
												std::unordered_set<int> &group_vertices, int &lambda)
{

	input_graph.clear();
	group_graph.clear();
	group_vertices.clear();

	std::string line_content;
	std::ifstream myfile(instance_name); // open the file
	if (myfile.is_open())				 // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			std::vector<std::string> Parsed_content = parse_string(line_content, " ");

			if (!Parsed_content[0].compare("input_graph") && !Parsed_content[1].compare("Vertex"))
			// when it's equal, compare returns 0
			{
				int v = std::stoi(Parsed_content[2]);
				int nw = std::stod(Parsed_content[3]);
				graph_hash_of_mixed_weighted_add_vertex(input_graph, v, nw);
			}
			else if (!Parsed_content[0].compare("input_graph") && !Parsed_content[1].compare("Edge"))
			{
				int v1 = std::stoi(Parsed_content[2]);
				int v2 = std::stoi(Parsed_content[3]);
				int ec = std::stod(Parsed_content[4]);
				graph_hash_of_mixed_weighted_add_edge(input_graph, v1, v2, ec);
			}
			else if (!Parsed_content[0].compare("group_graph") && !Parsed_content[1].compare("Vertex"))
			{
				int v = std::stoi(Parsed_content[2]);
				int nw = std::stod(Parsed_content[3]);
				graph_hash_of_mixed_weighted_add_vertex(group_graph, v, nw);
			}
			else if (!Parsed_content[0].compare("group_graph") && !Parsed_content[1].compare("Edge"))
			{
				int v1 = std::stoi(Parsed_content[2]);
				int v2 = std::stoi(Parsed_content[3]);
				int ec = std::stod(Parsed_content[4]);
				graph_hash_of_mixed_weighted_add_edge(group_graph, v1, v2, ec);
			}
			else if (!Parsed_content[0].compare("group_vertices"))
			{
				int g = std::stoi(Parsed_content[1]);
				group_vertices.insert(g);
			}
			else if (!Parsed_content[0].compare("lambda"))
			{
				lambda = std::stod(Parsed_content[1]);
			}
		}

		myfile.close(); // close the file
	}
	else
	{
		std::cout << "Unable to open file " << instance_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}
void graph_hash_of_mixed_weighted_read_community(std::string instance_name, int *community, int *c_size)
{

	std::string line_content;
	std::ifstream myfile(instance_name); // open the file
	getline(myfile, line_content);
	// std::vector<std::string> Parsed_content = parse_string(line_content, " ");
	//  for (size_t i = 0; i < 3; i++)
	//  {
	//  	c_size[i] = std::stod(Parsed_content[i]);
	//  }
	// printf("community size %d %d %d\n",c_size[0],c_size[1],c_size[2]);
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			std::vector<std::string> Parsed_content = parse_string(line_content, " ");
			int v = std::stoi(Parsed_content[0]);
			int c = std::stod(Parsed_content[1]);
			community[v] = c;
		}

		myfile.close(); // close the file
	}
	else
	{
		std::cout << "Unable to open file " << instance_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}

void graph_hash_of_mixed_weighted_read_for_GSTP_1(std::string instance_name,
												  graph_hash_of_mixed_weighted &input_graph, graph_hash_of_mixed_weighted &group_graph,
												  std::unordered_set<int> &group_vertices, int &lambda, std::vector<std::vector<int>> &community)
{

	input_graph.clear();
	group_graph.clear();
	group_vertices.clear();

	std::string line_content;
	std::ifstream myfile(instance_name); // open the file
	getline(myfile, line_content);
	std::vector<std::string> Parsed_content = parse_string(line_content, " ");
	int cnum = std::stod(Parsed_content[2]);
	cout << Parsed_content[0] << " " << Parsed_content[1] << endl;
	community.resize(3);
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			std::vector<std::string> Parsed_content = parse_string(line_content, " ");

			if (!Parsed_content[0].compare("input_graph") && !Parsed_content[1].compare("Vertex"))
			// when it's equal, compare returns 0
			{
				int v = std::stoi(Parsed_content[2]);
				int nw = std::stod(Parsed_content[3]);
				int c = std::stod(Parsed_content[4]);
				graph_hash_of_mixed_weighted_add_vertex(input_graph, v, nw);
				community[c].push_back(v);
			}
			else if (!Parsed_content[0].compare("input_graph") && !Parsed_content[1].compare("Edge"))
			{
				int v1 = std::stoi(Parsed_content[2]);
				int v2 = std::stoi(Parsed_content[3]);
				int ec = std::stod(Parsed_content[4]);
				graph_hash_of_mixed_weighted_add_edge(input_graph, v1, v2, ec);
			}
		}

		myfile.close(); // close the file
	}
	else
	{
		std::cout << "Unable to open file " << instance_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}

void graph_hash_of_mixed_weighted_read_for_Group(std::string instance_name,
												 graph_hash_of_mixed_weighted &input_graph, graph_hash_of_mixed_weighted &group_graph,
												 std::unordered_set<int> &group_vertices)
{

	std::string line_content;
	std::ifstream myfile(instance_name); // open the file
	int vnum = input_graph.hash_of_vectors.size();
	for (size_t i = 0; i < vnum; i++)
	{
		graph_hash_of_mixed_weighted_add_vertex(group_graph, i, 0);
	}
	
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			std::vector<std::string> Parsed_content = parse_string(line_content, ":");
			int g = std::stod(Parsed_content[0].substr(1, Parsed_content[0].length()))+vnum;
			graph_hash_of_mixed_weighted_add_vertex(group_graph, g, 0);
			std::vector<std::string> groups = parse_string(Parsed_content[1], " ");
			for (size_t i = 0; i < groups.size(); i++)
			{
				int v = std::stoi(groups[i]);
				graph_hash_of_mixed_weighted_add_edge(group_graph, g, v, 1);
			}
		}

		myfile.close(); // close the file
	}
	else
	{
		std::cout << "Unable to open file " << instance_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}

void graph_hash_of_mixed_weighted_read_for_inquire(std::string instance_name,
												   std::vector<std::vector<int>> &inquire)
{

	std::string line_content;
	std::ifstream myfile(instance_name); // open the file

	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			inquire.push_back({});
			std::vector<std::string> Parsed_content = parse_string(line_content, " ");
			for (size_t i = 0; i < Parsed_content.size()-1; i++)
			{	
				int v = std::stoi(Parsed_content[i]);	
				inquire[inquire.size()-1].push_back(v);
			}
		}

		myfile.close(); // close the file
	}
	else
	{
		std::cout << "Unable to open file " << instance_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}
