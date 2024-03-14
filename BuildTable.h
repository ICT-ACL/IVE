#pragma once
#include "graph.h"
#include "intersection_algos.hpp"
using namespace std;

struct BSRSet {
    int *base_ = nullptr;
    int *states_ = nullptr;
    int size_ = 0;

    BSRSet() = default;
    friend bool operator!=(BSRSet x, BSRSet y) {
        if (x.size_ != y.size_)
            return true;
        for (int i=0;i<x.size_;++i) {
            if (x.base_[i] != y.base_[i]) return true;
            if (x.states_[i] != y.states_[i]) return true;
        }
        return false;
    }
    friend bool operator<(BSRSet x, BSRSet y) {
        if (x.size_ != y.size_)
            return x.size_ < y.size_;
        for (int i = 0; i < x.size_; ++i) {
            if (x.base_[i] != y.base_[i]) return x.base_[i] < y.base_[i];
            if (x.states_[i] != y.states_[i]) return x.states_[i] < y.states_[i];
        }
        return false;
    }
};

struct BSRGraph {
    vector<BSRSet> bsrs;
    int max_d_ = 0;

    BSRGraph() = default;

    template<typename OFF, typename T>
    void load(size_t num_vertices, OFF &off_beg, OFF &node_off_end, T &adj) {
        bsrs.resize(num_vertices);

        int max_d = 0;
        {
            if (node_off_end != off_beg) {
                cout << "err" << endl;
                for (auto u = 0u; u < num_vertices; u++) {
                    node_off_end[u + 1] = static_cast<uint32_t>(
                            lower_bound(adj + off_beg[u], adj + off_beg[u + 1], u) - adj);
                }
            }
            for (auto u = 0u; u < num_vertices; u++) {
                max_d = max<int>(max_d, node_off_end[u + 1] - off_beg[u]);
            }
            int *tmp_base = new int[max_d];
            int *tmp_state = new int[max_d];

            for (int i = 0; i < num_vertices; i++) {

                auto degree = node_off_end[i + 1] - off_beg[i];
                auto tmp_size = offline_uint_trans_bsr(reinterpret_cast<int *>(adj) + off_beg[i], degree, tmp_base,
                                                       tmp_state);
                assert(tmp_size<= degree);
                bsrs[i].base_ = new int[degree + 16];
                bsrs[i].states_ = new int[degree + 16];
                bsrs[i].size_ = tmp_size;
                memcpy(bsrs[i].base_, tmp_base, static_cast<size_t>(tmp_size) * sizeof(int));
                memcpy(bsrs[i].states_, tmp_state, static_cast<size_t>(tmp_size) * sizeof(int));
            }

        }
        max_d_ = max_d;
    }
};

BSRGraph*** bsr_graph_;

void buildTables(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                            Edges ***edge_matrix, ui* order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui* flag = new ui[data_graph->getVerticesCount()];
    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    for (ui i = 0; i < query_vertices_num; ++i) {
        for (ui j = 0; j < query_vertices_num; ++j) 
            edge_matrix[i][j] = NULL;
    }

    std::vector<VertexID> build_table_order(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        build_table_order[i] = i;
    }

    std::sort(build_table_order.begin(), build_table_order.end(), [query_graph](VertexID l, VertexID r) {
        if (query_graph->getVertexDegree(l) == query_graph->getVertexDegree(r)) {
            return l < r;
        }
        return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
    });

    std::vector<ui> temp_edges(data_graph->getEdgesCount() * 2);

    for (auto u : build_table_order) {
        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

        ui updated_flag_count = 0;

        for (ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = u_nbrs[i];

            if (edge_matrix[u][u_nbr] != NULL)
                continue;

            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_nbr][u] = new Edges;
            edge_matrix[u_nbr][u]->vertex_count_ = candidates_count[u_nbr];
            edge_matrix[u_nbr][u]->offset_ = new ui[candidates_count[u_nbr] + 1];

            edge_matrix[u][u_nbr] = new Edges;
            edge_matrix[u][u_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_nbr]->offset_ = new ui[candidates_count[u] + 1];
            std::fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                VertexID v = candidates[u_nbr][j];
                edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count;
                const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

                ui local_degree = 0;

                for (ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = v_nbrs[k];

                    if (flag[v_nbr] != 0) {
                        ui position = flag[v_nbr] - 1;
                        temp_edges[local_edge_count++] = position;
                        edge_matrix[u][u_nbr]->offset_[position + 1] += 1;
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_nbr][u]->offset_[candidates_count[u_nbr]] = local_edge_count;
            edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_nbr][u]->edge_ = new ui[local_edge_count];
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

            edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_nbr]->edge_ = new ui[local_edge_count];

            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                }
                edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                ui begin = j;
                for (ui k = edge_matrix[u_nbr][u]->offset_[begin]; k < edge_matrix[u_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_nbr][u]->edge_[k];

                    edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end]++] = begin;
                }
            }

            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }
    
    int reverse_index[260];
    for (int i = 0; i < query_vertices_num; ++i) {
        int u = order[i];
        reverse_index[u] = i;
    }

    bsr_graph_ = new BSRGraph**[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        bsr_graph_[i] = new BSRGraph*[query_vertices_num];
        for (ui j = 0; j < query_vertices_num; ++j) {
            if (reverse_index[i] > reverse_index[j]) continue;

            bsr_graph_[i][j] = new BSRGraph[query_vertices_num];

            if (edge_matrix[i][j] != NULL) {
                bsr_graph_[i][j]->load(edge_matrix[i][j]->vertex_count_,
                                            edge_matrix[i][j]->offset_, edge_matrix[i][j]->offset_,
                                            edge_matrix[i][j]->edge_);
            }
        }
    }
}
