#pragma once
#include "pattern.h"
#include "prefix.h"

#include <cstdio>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <vector>

class Configuration {
public:
    Configuration(const std::string _adj_mat, int _size);
    ~Configuration();
    void init(Graph& complete);

    inline int get_total_prefix_num() const { return total_prefix_num; }
    inline int get_father_prefix_id(int prefix_id) const { return father_prefix_id[prefix_id]; }
    inline int get_loop_set_prefix_id(int loop) const { return loop_set_prefix_id[loop]; }
    inline int get_size() const { return size; }
    inline int get_last(int i) const { return last[i]; }
    inline int get_next(int i) const { return next[i]; }
    inline int get_in_exclusion_optimize_num() const { return in_exclusion_optimize_num; }
    inline void set_in_exclusion_optimize_num(int num) { in_exclusion_optimize_num = num; }
    inline int get_total_restrict_num() const { return total_restrict_num; }
    inline int get_restrict_last(int i) const { return restrict_last[i]; }
    inline int get_restrict_next(int i) const { return restrict_next[i]; }
    inline int get_restrict_index(int i) const { return restrict_index[i]; }
    inline long long get_in_exclusion_optimize_redundancy() const { return in_exclusion_optimize_redundancy; }
    inline const std::string get_adj_mat_ptr() const { return adj_mat; }
    inline void set_adj_mat(std::string _adj_mat) { adj_mat = _adj_mat; }

    void add_restrict(const std::vector<std::pair<int, int>> &restricts);
    void aggressive_optimize_get_all_pairs(std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector);
    void aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int>> isomorphism_vec, std::vector<std::vector<std::vector<int>>> permutation_groups, std::vector<std::pair<int, int>> ordered_pairs, std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector);

    void GraphZero_aggressive_optimize(std::vector<std::pair<int, int>> &ordered_pairs) const;
    void GraphZero_get_automorphisms(std::vector<std::vector<int>> &Aut) const;

    double our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    double GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt, unsigned int e_cnt);
    double Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &paris, int v_cnt, unsigned int e_cnt);

    void init_in_exclusion_optimize();

    std::vector<std::vector<std::vector<int>>> in_exclusion_optimize_group;
    std::vector<int> in_exclusion_optimize_val;
    std::vector<std::pair<int, int>> restrict_pair;

private:
    std::string adj_mat;
    std::vector<int> father_prefix_id;
    std::vector<int> last;
    std::vector<int> next;
    std::vector<int> loop_set_prefix_id;
    std::vector<Prefix> prefix;
    std::vector<int> restrict_last;
    std::vector<int> restrict_next;
    std::vector<int> restrict_index;
    int size;
    int total_prefix_num;
    int total_restrict_num;
    int in_exclusion_optimize_num;
    long long in_exclusion_optimize_redundancy;

    void build_loop_invariant();
    int find_father_prefix(int data_size, std::vector<int> data);

    void get_in_exclusion_optimize_group(int depth, std::vector<int> &id, int id_cnt, std::vector<int> &in_exclusion_val);
    //use principle of inclusion-exclusion to optimize
    void set_in_exclusion_optimize_redundancy(Graph& complete);
};
