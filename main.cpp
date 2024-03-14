#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cassert>
#include <string.h>
#include <bitset>
#include "matchingcommand.h"
#include "graph.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "reordering.h"
#include "enumeration.h"
using namespace std;

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_max_embedding_num = command.getMaximumEmbeddingNum();
    auto start = std::chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);

    Graph* data_graph = new Graph();
    data_graph->loadGraphFromFile(input_data_graph_file);

    ui** candidates = NULL;
    ui* candidates_count = NULL;
    Filter(data_graph, query_graph, candidates, candidates_count);
    sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    Edges ***edge_matrix = NULL;
    edge_matrix = new Edges **[query_graph->getVerticesCount()];
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
    }

    ui* matching_order = NULL;
    MDE(data_graph, query_graph, matching_order, candidates_count);
    buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix, matching_order);

    qfliter_bsr_graph_ = bsr_graph_;

    // print order
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        std::cout << matching_order[i] << " ";
    }
    std::cout << std::endl;

    size_t output_limit = 0;
    size_t embedding_count = 0;
    sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);

    size_t call_count = 0;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    auto start1= std::chrono::high_resolution_clock::now();
    embedding_count = explore(data_graph, query_graph, candidates, candidates_count,
        output_limit, call_count, matching_order);
    std::chrono::duration<double> end1 = std::chrono::high_resolution_clock::now() - start1;

    printf("#Embeddings: %zu\n", embedding_count);
    printf("Call Count: %zu\n", call_count);
    return 0;
}
