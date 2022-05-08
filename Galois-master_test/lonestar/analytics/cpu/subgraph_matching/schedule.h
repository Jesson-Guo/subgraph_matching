#pragma once
#include "pattern.h"
#include "prefix.h"
#include "disjoint_set_union.h"
#include "schedule.h"
#include "subgraph_matching.h"

#include <cstdio>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <vector>

class Schedule {
public:
    //TODO : more kinds of constructors to construct different Schedules from one Pattern
    Schedule(const Pattern& pattern, bool& is_pattern_valid, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize, Graph& complete, uint64_t max_degree, int v_cnt, unsigned int e_cnt, long long tri_cnt = 0);
    // performance_modeling type = 0 : not use modeling
    //                      type = 1 : use our modeling
    //                      type = 2 : use GraphZero's modeling
    //                      type = 3 : use naive modeling
    // restricts_type = 0 : not use restricts
    //                = 1 : use our restricts
    //                = 2 : use GraphZero's restricts
    Schedule(const int* _adj_mat, int _size, Graph& complete, uint64_t max_degree);
    ~Schedule();
    inline int get_total_prefix_num() const { return total_prefix_num; }
    inline int get_father_prefix_id(int prefix_id) const { return father_prefix_id[prefix_id]; }
    inline int get_loop_set_prefix_id(int loop) const { return loop_set_prefix_id[loop]; }
    inline int get_size() const { return size; }
    inline int get_last(int i) const { return last[i]; }
    inline int get_next(int i) const { return next[i]; }
    inline int get_in_exclusion_optimize_num() const { return in_exclusion_optimize_num; }
    inline void set_in_exclusion_optimize_num(int num) { in_exclusion_optimize_num = num; }
    void add_restrict(int* restrict_first, int* restrict_second, int restrict_size);
    inline int get_total_restrict_num() const { return total_restrict_num; }
    inline int get_restrict_last(int i) const { return restrict_last[i]; }
    inline int get_restrict_next(int i) const { return restrict_next[i]; }
    inline int get_restrict_index(int i) const { return restrict_index[i]; }
    inline long long get_in_exclusion_optimize_redundancy() const { return in_exclusion_optimize_redundancy; }
    int get_max_degree() const;
    int get_multiplicity() const;
    void aggressive_optimize(std::vector<std::pair<int, int>>& ordered_pairs) const;
    int aggressive_optimize_get_all_pairs(int** restricts_first, int** restricts_second, int* restricts_size);
    void aggressive_optimize_dfs(Pattern base_dag, int** isomorphism_vec, int iso_size, std::vector<std::vector<std::vector<int>>> permutation_groups, int* ordered_first, int* ordered_second, int ordered_size, int** restricts_first, int** restricts_second, int* restricts_size, int& rest_size);
    void restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt, std::vector<std::vector<std::pair<int, int>>> ordered_pairs_vector, std::vector<std::pair<int, int>> &best_restricts) const;
    int restricts_generate(const int* cur_adj_mat, int** restricts_first, int** restricts_second, int* restricts_size, Graph& complete, uint64_t max_degree);

    void GraphZero_aggressive_optimize(std::vector<std::pair<int, int>>& ordered_pairs) const;
    void GraphZero_get_automorphisms(std::vector<std::vector<int>>& Aut) const;

    int get_isomorphism_vec(int** isomorphism_vec) const;
    static std::vector<std::vector<int>> calc_permutation_group(int* vec, int size);
    inline const int* get_adj_mat_ptr() const { return adj_mat; }

    std::vector<std::vector<std::vector<int>>> in_exclusion_optimize_group;
    std::vector<int> in_exclusion_optimize_val;
    int* rest_first;
    int* rest_second;

    uint64_t max_degree;
    std::vector<std::vector<int>> isomorphism_vec;

private:
    int* adj_mat;
    int* father_prefix_id;
    int* last;
    int* next;
    int* loop_set_prefix_id;
    Prefix* prefix;
    int* restrict_last;
    int* restrict_next;
    int* restrict_index;
    int size;
    int total_prefix_num;
    int total_restrict_num;
    int in_exclusion_optimize_num;

    long long in_exclusion_optimize_redundancy;

    void build_loop_invariant();
    int find_father_prefix(int data_size, const int* data);
    void get_full_permutation(int** vec) const;

    double our_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    double GraphZero_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt);
    double Naive_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt);

    void get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val);
    //use principle of inclusion-exclusion to optimize
    void init_in_exclusion_optimize();

    int get_vec_optimize_num(int* vec);

    int remove_invalid_permutation(int** candidate_permutations, int n);

    void set_in_exclusion_optimize_redundancy(Graph& complete, uint64_t max_degree);
};
