#pragma once
#include "pattern.h"
#include "prefix.h"
#include "disjoint_set_union.h"

#include <assert.h>
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>
#include <atomic>
#include <queue>
#include <map>
#include <iostream>
#include <vector>

class Graph;
class Schedule;
class VertexSet;
class DataLoader;

class DataLoader {
public:

    bool load_data(Graph* &g, const char* path, int oriented_type = 0, long long tri_cnt=0);
        // pattern_diameter means max distance between two vertex in graph
        // oriented_type is used to reorder dataset
        // oriented_type == 0 do nothing
        //               == 1 high degree first
        //               == 2 low degree first
    bool load_complete(Graph* &g, int clique_size);

private:
    static bool cmp_pair(std::pair<int,int>a, std::pair<int,int>b);
    static bool cmp_degree_gt(std::pair<int,int> a,std::pair<int,int> b);
    static bool cmp_degree_lt(std::pair<int,int> a,std::pair<int,int> b);

    long long comb(int n,int k);
    bool general_load_data(Graph* &g, const char* path, int oriented_type = 0, long long tri_cnt = 0);

    std::map<int,int> id;
};

class VertexSet {
public:
    VertexSet();
    // allocate new memory according to max_intersection_size
    void init();
    void init(int init_size);
    // use memory from Graph, do not allocate new memory
    void init(int input_size, int* input_data);
    void copy(int input_size, const int* input_data);
    ~VertexSet();
    void intersection(const VertexSet& set0, const VertexSet& set1, int min_vertex = -1, bool clique = false);
    void intersection_with(const VertexSet& set1);
    //set1 is unordered
    static int unorderd_subtraction_size(const VertexSet& set0, const VertexSet& set1, int size_after_restrict = -1);
    void insert_ans_sort(int val);
    inline int get_size() const { return size;}
    inline int get_data(int i) const { return data[i];}
    inline const int* get_data_ptr() const { return data;}
    inline int* get_data_ptr() { return data;}
    inline void push_back(int val) { data[size++] = val;}
    inline void pop_back() { --size;}
    inline int get_last() const { return data[size - 1];}
    bool has_data(int val);
    static int max_intersection_size;
    void build_vertex_set(const Schedule& schedule, const VertexSet* vertex_set, int* input_data, int input_size, int prefix_id, int min_vertex = -1, bool clique = false);
private:
    int* data;
    int size;
    bool allocate;
};

class Schedule {
public:
    //TODO : more kinds of constructors to construct different Schedules from one Pattern
    Schedule(const Pattern& pattern, bool& is_pattern_valid, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize, int v_cnt, unsigned int e_cnt, long long tri_cnt = 0);
    // performance_modeling type = 0 : not use modeling
    //                      type = 1 : use our modeling
    //                      type = 2 : use GraphZero's modeling
    //                      type = 3 : use naive modeling
    // restricts_type = 0 : not use restricts
    //                = 1 : use our restricts
    //                = 2 : use GraphZero's restricts
    Schedule(const int* _adj_mat, int _size);
    ~Schedule();
    inline int get_total_prefix_num() const { return total_prefix_num;}
    inline int get_father_prefix_id(int prefix_id) const { return father_prefix_id[prefix_id];}
    inline int get_loop_set_prefix_id(int loop) const { return loop_set_prefix_id[loop];}
    inline int get_size() const { return size;}
    inline int get_last(int i) const { return last[i];}
    inline int get_next(int i) const { return next[i];}
    inline int get_in_exclusion_optimize_num() const { return in_exclusion_optimize_num;}
    inline void set_in_exclusion_optimize_num(int num) { in_exclusion_optimize_num = num; }
    int get_in_exclusion_optimize_num_when_not_optimize();
    void add_restrict(const std::vector< std::pair<int, int> >& restricts);
    inline int get_total_restrict_num() const { return total_restrict_num;}
    inline int get_restrict_last(int i) const { return restrict_last[i];}
    inline int get_restrict_next(int i) const { return restrict_next[i];}
    inline int get_restrict_index(int i) const { return restrict_index[i];}
    inline long long get_in_exclusion_optimize_redundancy() const { return in_exclusion_optimize_redundancy; }
    int get_max_degree() const;
    int get_multiplicity() const;
    void aggressive_optimize(std::vector< std::pair<int,int> >& ordered_pairs) const;
    void aggressive_optimize_get_all_pairs(std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector);
    void aggressive_optimize_dfs(Pattern base_dag, std::vector< std::vector<int> > isomorphism_vec, std::vector< std::vector< std::vector<int> > > permutation_groups, std::vector< std::pair<int,int> > ordered_pairs, std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector);
    void restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt, std::vector< std::vector< std::pair<int,int> > > ordered_pairs_vector, std::vector< std::pair<int,int> >& best_restricts) const;
    void restricts_generate(const int* cur_adj_mat, std::vector< std::vector< std::pair<int,int> > > &restricts);

    void GraphZero_aggressive_optimize(std::vector< std::pair<int,int> >& ordered_pairs) const;
    void GraphZero_get_automorphisms(std::vector< std::vector<int> > &Aut) const;

    std::vector< std::vector<int> > get_isomorphism_vec() const;
    static std::vector< std::vector<int> > calc_permutation_group(const std::vector<int> vec, int size);
    inline const int* get_adj_mat_ptr() const {return adj_mat;}
        

    void print_schedule() const;

    std::vector< std::vector< std::vector<int> > >in_exclusion_optimize_group;
    std::vector< int > in_exclusion_optimize_val;
    std::vector< std::pair<int,int> > restrict_pair;
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
    void get_full_permutation(std::vector< std::vector<int> >& vec, bool use[], std::vector<int> tmp_vec, int depth) const;
    void performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt);
    void bug_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt);
    void new_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    void GraphZero_performance_modeling(int* best_order, int v_cnt, unsigned int e_cnt);

    double our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    double GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt);
    double Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &paris, int v_cnt, unsigned int e_cnt);

    void get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val); 
    //use principle of inclusion-exclusion to optimize
    void init_in_exclusion_optimize();
    
    int get_vec_optimize_num(const std::vector<int> &vec);

    void remove_invalid_permutation(std::vector< std::vector<int> > &candidate_permutations);

    void set_in_exclusion_optimize_redundancy();
};

class Graph {
public:
    int v_cnt; // number of vertex
    unsigned int e_cnt; // number of edge
    long long tri_cnt; // number of triangle
    double max_running_time = 60 * 60 * 24; // second

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]
    
    Graph() {
        v_cnt = 0;
        e_cnt = 0;
        edge = nullptr;
        vertex = nullptr;
    }

    ~Graph() {
        if(edge != nullptr) delete[] edge;
        if(vertex != nullptr) delete[] vertex;
    }

    int intersection_size(int v1,int v2);
    int intersection_size_clique(int v1,int v2);

    //single thread triangle counting
    long long triangle_counting();
    
    //multi thread triangle counting
    long long triangle_counting_mt(int thread_count);

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule);

    //this function will be defined at code generation
    long long unfold_pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

    int max_degree;
private:
    void tc_mt(long long * global_ans);

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique = false);

    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    //this function will be defined at code generation
    void unfold_pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

};
