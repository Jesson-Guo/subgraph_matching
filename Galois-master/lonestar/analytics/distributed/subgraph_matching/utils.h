#pragma once
#include "vertex_set.h"

#include <assert.h>
#include <cstring>
#include <set>
#include <vector>
#include <cstdio>
#include <algorithm>

void get_full_permutation(int size, std::vector<std::vector<int>> &vec, bool use[], std::vector<int> tmp_vec, int depth) {
    if (depth == size) {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i) {
        if (use[i] == false) {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(size, vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
    }
}

// 计算置换群 like: (A)(B,C)(D)
std::vector<std::vector<int>> calc_permutation_group(const std::vector<int> vec, int size) {
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<std::vector<int>> res;
    res.clear();
    for (unsigned int i = 0; i < vec.size(); ++i) {
        if (use[i] == false) {
            std::vector<int> tmp_vec;
            tmp_vec.clear();
            tmp_vec.push_back(i);
            use[i] = true;
            int x = vec[i];
            while (use[x] == false) {
                use[x] = true;
                tmp_vec.push_back(x);
                x = vec[x];
            }
            res.push_back(tmp_vec);
        }
    }
    return res;
}

void aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int>> isomorphism_vec, std::vector<std::vector<std::vector<int>>> permutation_groups, std::vector<std::pair<int,int>> ordered_pairs, std::vector<std::vector<std::pair<int,int>>>& ordered_pairs_vector) {
    for (unsigned int i = 0; i < isomorphism_vec.size();) {
        Pattern test_dag(base_dag);
        const std::vector<int>& iso = isomorphism_vec[i];
        for (const std::pair<int, int>& pair : ordered_pairs)
            test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
        if (test_dag.is_dag() == false) {
            // is not dag means conflict
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
        }
        else ++i;
    }

    if (isomorphism_vec.size() == 1) {
        ordered_pairs_vector.push_back(ordered_pairs);
        return;
    }

    std::pair<int, int> found_pair;
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        for (const std::vector<int>& v : permutation_groups[i]) {
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int,int>(v[0], v[1]);
                std::vector<std::vector<std::vector<int>>> next_permutation_groups = permutation_groups;
                std::vector<std::vector<int>> next_isomorphism_vec = isomorphism_vec;
                std::vector<std::pair<int,int>> next_ordered_pairs = ordered_pairs;
                Pattern next_base_dag = base_dag;

                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                next_isomorphism_vec.erase(next_isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                next_ordered_pairs.push_back(found_pair);
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);

                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, next_permutation_groups, next_ordered_pairs, ordered_pairs_vector);
            }
        }
        if (two_element_number>= 1) break;
        else ++i;
    }
}

VertexSet intersection(const VertexSet &set0, const VertexSet &set1) {
    VertexSet set;
    if (&set0 == &set1) {
        return set0;
    }
    int i = 0;
    int j = 0;

    while (i < set0.get_size() && j < set1.get_size()) {
        data0 = set0.get_data(i);
        data1 = set1.get_data(j);
        if (data0 < data1) ++i;
        else if (data0 > data1) ++j;
        else {
            set.push_back(data0);
            ++i;
            ++j;
        }
    }
    return set;
}

int VertexSet::unorderd_subtraction_size(const VertexSet &set0, const VertexSet &set1, int size_after_restrict) {
    int size0 = set0.get_size();
    int size1 = set1.get_size();
    if (size_after_restrict != -1)
        size0 = size_after_restrict;

    int ret = size0;
    const std::vector<int> set_t = set0.get_data_ptr();
    for (int j = 0; j < size1; ++j)
        if (std::binary_search(set_t.begin(), set_t.begin() + size0, set1.get_data(j)))
            --ret;
    return ret;
}

void restricts_generate(const std::string cur_adj_mat, int size, std::vector<std::vector<std::pair<int,int>>> &restricts, Graph& complete, uint64_t max_degree) {
    Schedule schedule(cur_adj_mat, size, complete, max_degree);
    schedule.aggressive_optimize_get_all_pairs(restricts);
    long long ans = pattern_matching(complete, schedule, max_degree);
    ans /= schedule.get_multiplicity();
    for (int i = 0; i < (int)restricts.size();) {
        Schedule cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size(), complete, max_degree);
        cur_schedule.add_restrict(restricts[i]);
        long long cur_ans = pattern_matching(complete, cur_schedule, max_degree);
        if (cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        }
        else {
            ++i;
        }
    }
}
