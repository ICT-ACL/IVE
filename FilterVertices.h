#pragma once

#include "graph.h"
#include <vector>
#define INVALID_VERTEX_ID 1000000000


void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m) {
    int ptr;
    int i = 0;
    for(; i < n; i++) {
        int s_ptr = col_ptrs[i];
        int e_ptr = col_ptrs[i + 1];
        for(ptr = s_ptr; ptr < e_ptr; ptr++) {
            int r_id = col_ids[ptr];
            if(row_match[r_id] == -1) {
                match[i] = r_id;
                row_match[r_id] = i;
                break;
            }
        }
    }
}

void match_bfs(int* col_ptrs, int* col_ids, int* match, int* row_match, int* visited,
                        int* queue, int* previous, int n, int m) {
    int queue_ptr, queue_col, ptr, next_augment_no, i, j, queue_size,
            row, col, temp, eptr;

    old_cheap(col_ptrs, col_ids, match, row_match, n, m);

    memset(visited, 0, sizeof(int) * m);

    next_augment_no = 1;
    for(i = 0; i < n; i++) {
        if(match[i] == -1 && col_ptrs[i] != col_ptrs[i+1]) {
            queue[0] = i; queue_ptr = 0; queue_size = 1;

            while(queue_size > queue_ptr) {
                queue_col = queue[queue_ptr++];
                eptr = col_ptrs[queue_col + 1];
                for(ptr = col_ptrs[queue_col]; ptr < eptr; ptr++) {
                    row = col_ids[ptr];
                    temp = visited[row];

                    if(temp != next_augment_no && temp != -1) {
                        previous[row] = queue_col;
                        visited[row] = next_augment_no;

                        col = row_match[row];

                        if(col == -1) {
                            while(row != -1) {
                                col = previous[row];
                                temp = match[col];
                                match[col] = row;
                                row_match[row] = col;
                                row = temp;
                            }
                            next_augment_no++;
                            queue_size = 0;
                            break;
                        } else {
                            queue[queue_size++] = col;
                        }
                    }
                }
            }

            if(match[i] == -1) {
                for(j = 1; j < queue_size; j++) {
                    visited[match[queue[j]]] = -1;
                }
            }
        }
    }
}

bool verifyExactTwigIso(const Graph *data_graph, const Graph *query_graph, ui data_vertex, ui query_vertex,
                                   bool **valid_candidates, int *left_to_right_offset, int *left_to_right_edges,
                                   int *left_to_right_match, int *right_to_left_match, int* match_visited,
                                   int* match_queue, int* match_previous) {
    ui left_partition_size;
    ui right_partition_size;
    const VertexID* query_vertex_neighbors = query_graph->getVertexNeighbors(query_vertex, left_partition_size);
    const VertexID* data_vertex_neighbors = data_graph->getVertexNeighbors(data_vertex, right_partition_size);

    ui edge_count = 0;
    for (int i = 0; i < left_partition_size; ++i) {
        VertexID query_vertex_neighbor = query_vertex_neighbors[i];
        left_to_right_offset[i] = edge_count;

        for (int j = 0; j < right_partition_size; ++j) {
            VertexID data_vertex_neighbor = data_vertex_neighbors[j];

            if (valid_candidates[query_vertex_neighbor][data_vertex_neighbor]) {
                bool flag = true;
                left_to_right_edges[edge_count++] = j;
            }
        }
    }

    left_to_right_offset[left_partition_size] = edge_count;

    memset(left_to_right_match, -1, left_partition_size * sizeof(int));
    memset(right_to_left_match, -1, right_partition_size * sizeof(int));

    match_bfs(left_to_right_offset, left_to_right_edges, left_to_right_match, right_to_left_match,
                               match_visited, match_queue, match_previous, left_partition_size, right_partition_size);
    for (int i = 0; i < left_partition_size; ++i) {
        if (left_to_right_match[i] == -1)
            return false;
    }

    return true;
}

void computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                            ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
    const std::unordered_map<LabelID, ui>* query_vertex_nlf = query_graph->getVertexNLF(query_vertex);
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
    count = 0;

    for (ui j = 0; j < data_vertex_num; ++j) {
        ui data_vertex = data_vertices[j];
        if (data_graph->getVertexDegree(data_vertex) >= degree) {

            // NFL check
            const std::unordered_map<LabelID, ui>* data_vertex_nlf = data_graph->getVertexNLF(data_vertex);

            if (data_vertex_nlf->size() >= query_vertex_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_nlf) {
                    auto iter = data_vertex_nlf->find(element.first);
                    if (iter == data_vertex_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    if (buffer != NULL) {
                        buffer[count] = data_vertex;
                    }
                    count += 1;
                }
            }
        }
    }
}

void sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}

void allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                ui *&candidates_count) {
    ui query_vertex_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[query_vertex_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui*[query_vertex_num];

    for (ui i = 0; i < query_vertex_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}
void compactCandidates(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != INVALID_VERTEX_ID) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
    }
}
bool isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}


bool Filter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        computeCandidateWithNLF(data_graph, query_graph, query_vertex, candidates_count[query_vertex], candidates[query_vertex]);
    }

    // Allocate buffer.
    ui query_vertex_num = query_graph->getVerticesCount();
    ui data_vertex_num = data_graph->getVerticesCount();

    bool** valid_candidates = new bool*[query_vertex_num];
    for (ui i = 0; i < query_vertex_num; ++i) {
        valid_candidates[i] = new bool[data_vertex_num];
        memset(valid_candidates[i], 0, sizeof(bool) * data_vertex_num);
    }

    ui query_graph_max_degree = query_graph->getGraphMaxDegree();
    ui data_graph_max_degree = data_graph->getGraphMaxDegree();

    int* left_to_right_offset = new int[query_graph_max_degree*2 + 1];
    int* left_to_right_edges = new int[query_graph_max_degree*2 * data_graph_max_degree*2];
    int* left_to_right_match = new int[query_graph_max_degree*2];
    int* right_to_left_match = new int[data_graph_max_degree*2];
    int* match_visited = new int[data_graph_max_degree*2 + 1];
    int* match_queue = new int[query_vertex_num*2];
    int* match_previous = new int[data_graph_max_degree*2 + 1];

    int*** last_matches = new int**[query_vertex_num];
    // Record valid candidate vertices for each query vertex.
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        last_matches[i] = new int*[candidates_count[i]];
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];
            valid_candidates[query_vertex][data_vertex] = true;
        }
    }

    // Global refinement.
    int h=0, t=0;
    int *queue = new int[query_vertex_num+1];
    bool *change_set = new bool[query_vertex_num];
    for(int i=0;i<query_vertex_num;++i){
        change_set[i] = true;
        queue[h++] = i;
    }
    while(h!=t) {
        bool first = false;
        VertexID query_vertex = queue[t++];
        if (t == query_vertex_num + 1)
            t = 0;
        change_set[query_vertex] = false;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];
            if (data_vertex == INVALID_VERTEX_ID)
                continue;
            if (!verifyExactTwigIso(data_graph, query_graph, data_vertex, query_vertex, valid_candidates,
                                    left_to_right_offset, left_to_right_edges, left_to_right_match,
                                    right_to_left_match, match_visited, match_queue, match_previous)) {
                candidates[query_vertex][j] = INVALID_VERTEX_ID;
                valid_candidates[query_vertex][data_vertex] = false;
                ui u1_neighbors;
                const VertexID* query_vertex_neighbors = query_graph->getVertexNeighbors(query_vertex, u1_neighbors);
                if (!first){
                    for (int k=0;k<u1_neighbors;++k) if (!change_set[query_vertex_neighbors[k]]) {
                        change_set[query_vertex_neighbors[k]] = true;
                        queue[h++] = query_vertex_neighbors[k];
                        if (h==query_vertex_num+1) 
                            h = 0;
                    }
                    first = true;
                }
            }
        }
    }
    // Compact candidates.
    compactCandidates(candidates, candidates_count, query_vertex_num);

    // Release memory.
    for (ui i = 0; i < query_vertex_num; ++i) {
        delete[] valid_candidates[i];
    }

    delete[] valid_candidates;
    delete[] left_to_right_offset;
    delete[] left_to_right_edges;
    delete[] left_to_right_match;
    delete[] right_to_left_match;
    delete[] match_visited;
    delete[] match_queue;
    delete[] match_previous;

    return isCandidateSetValid(candidates, candidates_count, query_vertex_num);
}
