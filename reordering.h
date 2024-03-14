#pragma once
#include "graph.h"
#include <vector>

void MDE(const Graph *data_graph, const Graph *query_graph, ui *&order, ui *candidates_count) {
    bool is_isolated[260];
    int query_vertex_number = query_graph->getVerticesCount();
    order = new ui[query_vertex_number*2];
    bool used[260];
    bool extendable[260];
    for (int i = 0; i < query_vertex_number; ++i) {
        extendable[i] = false;
        is_isolated[i] = false;
        used[i] = false;
    }

    int start_vertex = 0;
    for (int i = 1; i < query_vertex_number; ++i) {
        if (query_graph->getVertexDegree(i) > query_graph->getVertexDegree(start_vertex))
            start_vertex = i;
        else if (query_graph->getVertexDegree(i) == query_graph->getVertexDegree(start_vertex)) {
            if (candidates_count[i] < candidates_count[start_vertex])
                start_vertex = i;
        }
    }
    order[0] = start_vertex;
    used[start_vertex] = true;

    ui edge_count;
    const ui* neighbors = query_graph -> getVertexNeighbors(start_vertex, edge_count);
    for (int i = 0; i < edge_count; ++i)
        extendable[neighbors[i]] = true;
    for (int loop = 1; loop < query_vertex_number; ++loop) {
        int best = -1;
        int scoreA;
        for (int i = 0; i < query_vertex_number; ++i) if (extendable[i]) {
            int ScoreB=0;
            ui edge_count;
            const ui* neighbors = query_graph -> getVertexNeighbors(i, edge_count);
            for (int p = 0; p < edge_count; ++p) {
                int un = neighbors[p];
                if (!used[un])
                    ++ScoreB;
            }
            if (ScoreB == 0) {
                best = i;
                scoreA = ScoreB;
                break;
            }

            if (best == -1) {
                best = i;
                scoreA = ScoreB;
                continue;
            }
            if (ScoreB > scoreA || (ScoreB == scoreA && candidates_count[i] < candidates_count[best])) {
                best = i;
                scoreA = ScoreB;
            }
        }
        if (scoreA == 0) 
            is_isolated[best] = true;
        assert(extendable[best]);
        extendable[best] = false;
        order[loop] = best;
        used[best] = true;
        ui edge_count;
        const ui* neighbors = query_graph -> getVertexNeighbors(best, edge_count);
        for (int i = 0; i < edge_count; ++i)
            if (!used[neighbors[i]])
                extendable[neighbors[i]] = true;
    }

    vector<int> tmp(query_vertex_number);
    int l = 0, r = query_vertex_number;
    for (int i = 0; i < query_vertex_number;++i) {
        int v = order[i];
        if (!is_isolated[v])
            tmp[l++] = v;
        else tmp[--r] = v;
    }
    assert(l == r);
    for (int i = 0; i < query_vertex_number;++i) 
        order[i] = tmp[i];
    printf("iso: %d\n", r);
}
