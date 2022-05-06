#pragma once
#include "vertex_set.h"

#include <assert.h>
#include <cstring>
#include <set>
#include <vector>
#include <cstdio>
#include <algorithm>

#ifndef INDEX
#define INDEX(x, y, n) ((x) * (n) + (y))
#endif

int get_max_degree(int size, std::string adj_mat) {
    int mx = 0;
    for (int i = 0; i < size; ++i) {
        int cnt = 0;
        for (int j = 0; j < size; ++j) cnt += adj_mat[INDEX(i,j,size)] - '0';
        if (cnt > mx) mx = cnt;
    }
    return mx;
}


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

// the i's vertex must connect to at lease one of i-1 vertices at front
void remove_invalid_permutation(std::vector<std::vector<int>>& candidate_permutations, int size, std::string adj_mat) {
    for (unsigned int i = 0; i < candidate_permutations.size();) {
        const auto& vec = candidate_permutations[i];
        bool tag = true;
        for (int x = 1; x < size; ++x) {
            bool have_edge = false;
            for (int y = 0; y < x; ++y) {
                if (adj_mat[INDEX(vec[x],vec[y],size)] == '1') {
                    have_edge = true;
                    break;
                }
            }
            if (!have_edge) {
                tag = false;
                break;
            }
        }
        if (tag) {
            ++i;
        }
        else {
            candidate_permutations.erase(candidate_permutations.begin() + i);
        }
    }
}

std::vector<std::vector<int>> get_isomorphism_vec(int size, std::string adj_mat) {
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    std::vector<std::vector<int>> vec;
    vec.clear();
    bool use[size];
    for (int i = 0; i < size; ++i) use[i] = false;
    std::vector<int> tmp_vec;
    get_full_permutation(size, vec, use, tmp_vec, 0);
    assert(vec.size() == pow);
    std::vector<std::vector<int>> isomorphism_vec;
    isomorphism_vec.clear();
    for (const std::vector<int>& v : vec) {
        bool flag = true;
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (adj_mat[INDEX(i, j, size)] != '0') {
                    if (adj_mat[INDEX(v[i], v[j], size)] == '0') {
                        // not isomorphism
                        flag = false;
                        break;
                    }
                }
            }
        }
        if (flag == true) isomorphism_vec.push_back(v);
    }
    return isomorphism_vec;
}

int get_multiplicity(int size, std::string adj_mat) {
    std::vector<std::vector<int>> isomorphism_vec = get_isomorphism_vec(size, adj_mat);
    return isomorphism_vec.size();
}

int get_vec_optimize_num(int size, const std::vector<int>& vec, std::string adj_mat) {
    bool is_valid = true;
    for (int i = 1; i < size; ++i) {
        bool have_edge = false;
        for (int j = 0; j < i; ++j) {
            if (adj_mat[INDEX(vec[i], vec[j], size)] == '1') {
                have_edge = true;
                break;
            }
        }
        if (have_edge == false) {
            is_valid = false;
            break;
        }
    }
    if (!is_valid) return -1;

    for (int k = 2; k <= size; ++k) {
        bool flag = true;
        for (int i = size - k + 1; i < size; ++i) {
            if (adj_mat[INDEX(vec[size - k], vec[i], size)] == '1') {
                flag = false;
                break;
            }
        }
        if (flag == false) return k - 1;
    }
    assert(0);
    return -1;
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

//set1 is unordered
int unorderd_subtraction_size(const VertexSet &set0, const VertexSet &set1, int size_after_restrict=-1) {
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

void print_schedule(int size, std::string adj_mat) {
    printf("Schedule:\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) printf("%c", adj_mat[INDEX(i,j,size)]);
        puts("");
    }
}
