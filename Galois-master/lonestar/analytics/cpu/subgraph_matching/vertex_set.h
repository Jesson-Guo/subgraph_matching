#pragma once
#include <algorithm>

class VertexSet {
public:
    VertexSet();
    ~VertexSet();
    // allocate new memory according to max_intersection_size
    void init();
    void init(int init_size);
    // use memory from Graph, do not allocate new memory
    void init(int input_size, std::vector<int> input_data);

    inline int get_size() const { return size; }
    inline int get_data(int i) const { return data[i]; }
    inline const std::vector<int> get_data_ptr() const { return data; }
    inline std::vector<int> get_data_ptr() { return data; }
    inline void push_back(int val) { data[size++] = val; }
    inline void pop_back() { --size; }
    inline int get_last() const { return data[size - 1]; }

    void intersection(const VertexSet set0, const VertexSet set1);
    bool has_data(int val);
    static int max_intersection_size;
    void build_vertex_set(int father_id, const std::vector<VertexSet> vertex_set, std::vector<int> input_data, int input_size);

private:
    std::vector<int> data;
    int size;
    bool allocate;
};

int VertexSet::max_intersection_size = -1;

VertexSet::VertexSet() : data(std::vector<int>(0)), size(0), allocate(false) {}

void VertexSet::init() {
    if (allocate == true)
        size = 0; // do not reallocate
    else {
        size = 0;
        allocate = true;
        data = std::vector<int>(max_intersection_size);
    }
}

//this function is only used for tmp_set in graph.cpp (i.e., init(Graph.max_degree))
void VertexSet::init(int input_size) {
    if (allocate == true)
        size = 0;
    else {
        size = 0;
        allocate = true;
        data = std::vector<int>(input_size);
    }
}

void VertexSet::init(int input_size, std::vector<int> input_data) {
    size = input_size;
    data = input_data;
    allocate = false;
}

VertexSet::~VertexSet() {}

void VertexSet::intersection(const VertexSet set0, const VertexSet set1) {
    if (&set0 == &set1) {
        size = set0.get_size();
        data = set0.get_data_ptr();
        return;
    }
    int i = 0;
    int j = 0;
    while (i < set0.get_size() && j < set1.get_size()) {
        int data0 = set0.get_data(i);
        int data1 = set1.get_data(j);
        if (data0 < data1) ++i;
        else if (data0 > data1) ++j;
        else {
            push_back(data0);
            ++i;
            ++j;
        }
    }
}

void VertexSet::build_vertex_set(int father_id, const std::vector<VertexSet> vertex_set, std::vector<int> input_data, int input_size) {
    // int father_id = schedule.get_father_prefix_id(prefix_id);
    if (father_id == -1)
        init(input_size, input_data);
    else {
        init();
        VertexSet tmp_vset;
        tmp_vset.init(input_size, input_data);
        intersection(vertex_set[father_id], tmp_vset);
    }
}

bool VertexSet::has_data(int val) {
    for (int i = 0; i < size; ++i)
        if (data[i] == val)
            return true;
    return false;
}
