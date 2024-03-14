#pragma once

#include<bits/stdc++.h>

typedef unsigned int ui;

typedef uint32_t VertexID;
typedef ui LabelID;

class Edges {
public:
    ui* offset_;
    ui* edge_;
    ui vertex_count_;
    ui edge_count_;
    ui max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

class Graph {
private:
    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;

    ui* offsets_;
    VertexID * neighbors_;
    LabelID* labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    std::unordered_map<LabelID, ui> labels_frequency_;

    ui* labels_offsets_;
    std::unordered_map<LabelID, ui>* nlf_;

private:
    void BuildReverseIndex() {
        reverse_index_ = new ui[vertices_count_];
        reverse_index_offsets_= new ui[labels_count_ + 1];
        reverse_index_offsets_[0] = 0;

        ui total = 0;
        for (ui i = 0; i < labels_count_; ++i) {
            reverse_index_offsets_[i + 1] = total;
            total += labels_frequency_[i];
        }

        for (ui i = 0; i < vertices_count_; ++i) {
            LabelID label = labels_[i];
            reverse_index_[reverse_index_offsets_[label + 1]++] = i;
        }
    }
    void BuildNLF() {
        nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];
        for (ui i = 0; i < vertices_count_; ++i) {
            ui count;
            const VertexID * neighbors = getVertexNeighbors(i, count);

            for (ui j = 0; j < count; ++j) {
                VertexID u = neighbors[j];
                LabelID label = getVertexLabel(u);
                if (nlf_[i].find(label) == nlf_[i].end()) {
                    nlf_[i][label] = 0;
                }

                nlf_[i][label] += 1;
            }
        }
    }

public:
    Graph() {
        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        labels_frequency_.clear();

        labels_offsets_ = NULL;
        nlf_ = NULL;
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] labels_offsets_;
        delete[] nlf_;
    }
    void loadGraphFromFile(const std::string &file_path) {
        std::ifstream infile(file_path);

        if (!infile.is_open()) {
            std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
            exit(-1);
        }

        char type;
        infile >> type >> vertices_count_ >> edges_count_;
        offsets_ = new ui[vertices_count_ +  1];
        offsets_[0] = 0;

        neighbors_ = new VertexID[edges_count_ * 2];
        labels_ = new LabelID[vertices_count_];
        labels_count_ = 0;
        max_degree_ = 0;

        LabelID max_label_id = 0;
        std::vector<ui> neighbors_offset(vertices_count_, 0);

        while (infile >> type) {
            if (type == 'v') { // Read vertex.
                VertexID id;
                LabelID  label;
                ui degree;
                infile >> id >> label >> degree;

                labels_[id] = label;
                offsets_[id + 1] = offsets_[id] + degree;

                if (degree > max_degree_) {
                    max_degree_ = degree;
                }

                if (labels_frequency_.find(label) == labels_frequency_.end()) {
                    labels_frequency_[label] = 0;
                    if (label > max_label_id)
                        max_label_id = label;
                }

                labels_frequency_[label] += 1;
            }
            else if (type == 'e') { // Read edge.
                VertexID begin;
                VertexID end;
                infile >> begin >> end;

                ui offset = offsets_[begin] + neighbors_offset[begin];
                neighbors_[offset] = end;

                offset = offsets_[end] + neighbors_offset[end];
                neighbors_[offset] = begin;

                neighbors_offset[begin] += 1;
                neighbors_offset[end] += 1;
            }
        }

        infile.close();
        labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

        for (auto element : labels_frequency_) {
            if (element.second > max_label_frequency_) {
                max_label_frequency_ = element.second;
            }
        }

        for (ui i = 0; i < vertices_count_; ++i) {
            std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
        }

        BuildReverseIndex();
        BuildNLF();
    }

    void printGraphMetaData() {
        std::cout << "|V|: " << vertices_count_ << ", |E|: " << edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
        std::cout << "Max Degree: " << max_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
    }

    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }


    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const ui * getNeighborsByLabel(const VertexID id, const LabelID label, ui& count) const {
        ui offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return nlf_ + id;
    }

    bool checkEdgeExistence(const VertexID u, const VertexID v, const LabelID u_label) const {
        ui count = 0;
        const VertexID* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
};
