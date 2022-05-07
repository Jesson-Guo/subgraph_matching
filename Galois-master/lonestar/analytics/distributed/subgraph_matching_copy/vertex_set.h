#pragma once
#include "schedule.h"
#include <algorithm>

class VertexSet {
public:
    VertexSet();
    // allocate new memory according to max_intersection_size
    void init();
    void init(int init_size);
    // use memory from Graph, do not allocate new memory
    void init(int input_size, std::vector<int> input_data);
    void copy(int input_size, const std::vector<int> input_data);
    ~VertexSet();
    void intersection(const VertexSet &set0, const VertexSet &set1, int min_vertex = -1, bool clique = false);
    void intersection_with(const VertexSet &set1);
    //set1 is unordered
    static int unorderd_subtraction_size(const VertexSet &set0, const VertexSet &set1, int size_after_restrict = -1);
    void insert_ans_sort(int val);
    inline int get_size() const { return size; }
    inline int get_data(int i) const { return data[i]; }
    inline const std::vector<int> get_data_ptr() const { return data; }
    inline std::vector<int> get_data_ptr() { return data; }
    inline void push_back(int val) { data[size++] = val; }
    inline void pop_back() { --size; }
    inline int get_last() const { return data[size - 1]; }
    bool has_data(int val);
    static int max_intersection_size;
    void build_vertex_set(const Schedule &schedule, const std::vector<VertexSet> vertex_set, std::vector<int> input_data, int input_size, int prefix_id, int min_vertex = -1, bool clique = false);

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

VertexSet::~VertexSet() {
    if (allocate== true) data.clear();
}

void VertexSet::intersection(const VertexSet &set0, const VertexSet &set1, int min_vertex, bool clique) {
    if (&set0 == &set1) {
        size = set0.get_size();
        data = set0.get_data_ptr();
        return;
    }
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    if (clique == true)
        if (set0.get_data(0) >= min_vertex || set1.get_data(0) >= min_vertex)
            return;
    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    if (clique) {
        // TODO : Try more kinds of calculation.
        // For example, we can use binary search find the last element which is smaller than min_vertex, and set its index as loop_size.
        while (i < size0 && j < size1) {
            if (data0 < data1) {
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
            }
            else if (data0 > data1) {
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            }
            else {
                push_back(data0);
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            }
        }
    }
    else {
        while (i < size0 && j < size1) {
            data0 = set0.get_data(i);
            data1 = set1.get_data(j);
            if (data0 < data1)
                ++i;
            else if (data0 > data1)
                ++j;
            else {
                push_back(data0);
                ++i;
                ++j;
            }
        }
    }
}

void VertexSet::intersection_with(const VertexSet &set1) {
    if (this == &set1)
        return;
    const VertexSet &set0 = *this;
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    size = 0;
    while (i < size0 && j < size1) {
        data0 = set0.get_data(i);
        data1 = set1.get_data(j);
        if (data0 < data1)
            ++i;
        else if (data0 > data1)
            ++j;
        else {
            push_back(data0);
            ++i;
            ++j;
        }
    }
}

void VertexSet::build_vertex_set(const Schedule &schedule, const std::vector<VertexSet> vertex_set, std::vector<int> input_data, int input_size, int prefix_id, int min_vertex, bool clique) {
    int father_id = schedule.get_father_prefix_id(prefix_id);
    if (father_id == -1)
        init(input_size, input_data);
    else {
        init();
        VertexSet tmp_vset;
        tmp_vset.init(input_size, input_data);
        intersection(vertex_set[father_id], tmp_vset, min_vertex, clique);
    }
}

bool VertexSet::has_data(int val) {
    for (int i = 0; i < size; ++i)
        if (data[i] == val)
            return true;
    return false;
}

int VertexSet::unorderd_subtraction_size(const VertexSet &set0, const VertexSet &set1, int size_after_restrict) {
    int size0 = set0.get_size();
    int size1 = set1.get_size();
    if (size_after_restrict != -1)
        size0 = size_after_restrict;

    int ret = size0;
    const std::vector<int> set0_ptr = set0.get_data_ptr();
    for (int j = 0; j < size1; ++j)
        if (std::binary_search(set0_ptr.begin(), set0_ptr.begin() + size0, set1.get_data(j)))
            --ret;
    return ret;
}
