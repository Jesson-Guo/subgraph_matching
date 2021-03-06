#pragma once
#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "Lonestar/BoilerPlate.h"
#include "Lonestar/Utils.h"

struct NodeData {
    int degree;
    std::vector<int> adj_node;
};

constexpr static const unsigned CHUNK_SIZE = 128U;

//! [withnumaalloc]
using Graph = galois::graphs::LC_CSR_Graph<NodeData, void>::
    with_no_lockable<true>::type ::with_numa_alloc<true>::type;
//! [withnumaalloc]
typedef Graph::GraphNode GNode;

#include "prefix.h"
#include "pattern.h"
#include "schedule.h"
#include "vertex_set.h"

#include <iostream>
#include <algorithm>

void pattern_matching_aggressive_func(Graph& graph, const Schedule &schedule, VertexSet* vertex_set, VertexSet &subtraction_set, VertexSet &tmp_set, uint64_t max_degree, galois::GAccumulator<long long>& local_ans, int depth) {
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth); // @@@
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;

    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    //Case: in_exclusion_optimize_num > 1
    if (depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num()) {
        int in_exclusion_optimize_num = schedule.get_in_exclusion_optimize_num(); // @@@
        int loop_set_prefix_ids[in_exclusion_optimize_num];
        loop_set_prefix_ids[0] = loop_set_prefix_id;
        for (int i = 1; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id(depth + i);
        for (int optimize_rank = 0; optimize_rank < (int)schedule.in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector<std::vector<int>> &cur_graph = schedule.in_exclusion_optimize_group[optimize_rank];
            long long val = schedule.in_exclusion_optimize_val[optimize_rank];
            for (int cur_graph_rank = 0; cur_graph_rank < (int)cur_graph.size(); ++cur_graph_rank) {
                //if size == 1 , we will not call intersection(...)
                //so we will not allocate memory for data
                //otherwise, we need to copy the data to do intersection(...)
                if (cur_graph[cur_graph_rank].size() == 1) {
                    int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    val = val * VertexSet::unorderd_subtraction_size(vertex_set[id], subtraction_set);
                }
                else {
                    int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
                    tmp_set.init(max_degree);
                    tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

                    for (int i = 2; i < (int)cur_graph[cur_graph_rank].size(); ++i) {
                        int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
                        tmp_set.intersection_with(vertex_set[id]);
                    }
                    val = val * VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                }
                if (val == 0)
                    break;
            }
            local_ans += val;
        }
        return;
    }
    //Case: in_exclusion_optimize_num <= 1
    if (depth == schedule.get_size() - 1) {
        // TODO : try more kinds of calculation. @@@
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0) {
            int min_vertex = graph.size();
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            const VertexSet &vset = vertex_set[loop_set_prefix_id];
            int size_after_restrict = std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex) - vset.get_data_ptr();
            if (size_after_restrict > 0) {
                // ???????????????????????????????????????
                local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
            }
        }
        else
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    // TODO : min_vertex is also a loop invariant
    int min_vertex = graph.size();
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
    int ii = 0;
    for (int &i = ii; i < loop_size; ++i) {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id, vertex);
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if (is_zero)
            continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        pattern_matching_aggressive_func(graph, schedule, vertex_set, subtraction_set, tmp_set, max_degree, local_ans, depth + 1);
        subtraction_set.pop_back();
    }
}

long long pattern_matching(Graph& graph, const Schedule &schedule, uint64_t max_degree) {
    VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
    VertexSet subtraction_set;
    VertexSet tmp_set;
    subtraction_set.init();
    galois::GAccumulator<long long> local_ans;
    // TODO : try different chunksize
    galois::do_all(
        galois::iterate(graph),
        [&](uint64_t vertex) {
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(
                    schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id);
            }
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(graph, schedule, vertex_set, subtraction_set, tmp_set, max_degree, local_ans, 1);
            subtraction_set.pop_back();
        },
        galois::loopname("pattern matching"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    delete[] vertex_set;
    std::cout << "local ans: " << local_ans.reduce() << std::endl;
    long long global_ans = local_ans.reduce();
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

// Schedule functions
Schedule::Schedule(const Pattern& pattern, bool &is_pattern_valid, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize, Graph& complete, uint64_t _max_degree, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    if( performance_modeling_type != 0 && tri_cnt == -1) {
        printf("Fatal: Can not use performance modeling if not have triangle number of this dataset.\n");
        fflush(stdout);
        assert(0);
    }

    is_pattern_valid = true;
    size = pattern.get_size();
    adj_mat = new int[size * size];
    max_degree = _max_degree;

    // not use performance_modeling, simply copy the adj_mat from pattern
    memcpy(adj_mat, pattern.get_adj_mat_ptr(), size * size * sizeof(int));

    std::vector< std::pair<int,int> > best_pairs;
    best_pairs.clear();
    //Initialize adj_mat
    //If we use performance_modeling, we may change the order of vertex,
    //the best order produced by performance_modeling(...) is saved in best_order[]
    //Finally, we use best_order[] to relocate adj_mat
    if (performance_modeling_type != 0) { 
        int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;

        int** candidate_permutations_t = new int*[pow];
        int** candidate_permutations = new int*[pow];
        for(int i=0; i<pow; i++) candidate_permutations[i] = new int[size];

        get_full_permutation(candidate_permutations);

        int perm_size = remove_invalid_permutation(candidate_permutations, pow);

        if(performance_modeling_type == 1) {
            //reduce candidates
            int max_val = 0;
            for(int i=0; i<perm_size; i++) {
                int* vec = new int[size];
                memcpy(vec, candidate_permutations[i], sizeof(int)*size);
                max_val = std::max(max_val, get_vec_optimize_num(vec));
                delete[] vec;
            }
            int j = 0;
            for(int i=0; i<perm_size; i++) {
                int* vec = new int[size];
                memcpy(vec, candidate_permutations[i], sizeof(int)*size);
                if (get_vec_optimize_num(vec) == max_val) {
                    candidate_permutations_t[j] = new int[size];
                    memcpy(candidate_permutations_t[j], vec, sizeof(int)*size);
                    j++;
                }
                delete[] vec;
            }
            perm_size = j;
        }
        for(int i=0;i<pow;i++) delete[] candidate_permutations[i];
        delete[] candidate_permutations;

        int *best_order = new int[size];
        double min_val = 0;
        bool have_best = false;

        // galois::do_all(
        //     galois::iterate(candidate_permutations),
        //     [&](auto &vec) {
        //     },
        //     galois::loopname("max intersection size calculate"),
        //     galois::chunk_size<CHUNK_SIZE>(),
        //     galois::steal(),
        //     galois::no_stats()
        // );

        // for(const auto &vec : candidate_permutations) {
        for(int o=0; o<perm_size; o++) {
            int* vec = new int[size];
            memcpy(vec, candidate_permutations_t[o], sizeof(int)*size);
            int* rank = new int[size];
            // int rank[size];
            for(int i = 0; i < size; ++i) rank[vec[i]] = i;

            int* cur_adj_mat = new int[size*size];
            for(int i = 0; i < size; ++i)
                for(int j = 0; j < size; ++j)
                    cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

            int pow = 1;
            for (int i = 2; i <= size; ++i) pow *= i;
            int** restricts_first = new int*[pow];
            int** restricts_second = new int*[pow];
            for (int j = 0; j < size-1; j++) {
                restricts_first[j] = new int[size-1];
                restricts_second[j] = new int[size-1];
            }
            int* restricts_size = new int[pow];
            int rest_size = 0;

            if(restricts_type == 1) {
                rest_size = restricts_generate(cur_adj_mat, restricts_first, restricts_second, restricts_size, complete, max_degree);
            }

            if(rest_size == 0) {
                std::vector< std::pair<int,int> > Empty;
                Empty.clear();

                double val;
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                    else {
                        val = Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                }

                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }

            // for(const auto& pairs : restricts_vector) {
            for(int i=0; i<rest_size; i++) {
                double val;
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, restricts_first[i], restricts_second[i], v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, restricts_first[i], restricts_second[i], v_cnt, e_cnt);
                    }
                    else {
                        val = Naive_estimate_schedule_restrict(vec, restricts_first[i], restricts_second[i], v_cnt, e_cnt);
                    }
                }

                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }
            delete[] rank;
            delete[] cur_adj_mat;
        }

        int rank[size];
        for(int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const int* pattern_adj_mat = pattern.get_adj_mat_ptr();
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)]; 
        delete[] best_order;
    }
    else {
        int* I = new int[size];
        for(int i = 0; i < size; ++i) I[i] = i;

        int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;
        int** restricts_first = new int*[pow];
        int** restricts_second = new int*[pow];
        for (int j = 0; j < size-1; j++) {
            restricts_first[j] = new int[size-1];
            restricts_second[j] = new int[size-1];
        }
        int* restricts_size = new int[pow];
        int rest_size = 0;

        if(restricts_type != 0) {
            if(restricts_type == 1) {
                rest_size = restricts_generate(adj_mat, restricts_first, restricts_second, restricts_size, complete, max_degree);
            }
        }

        bool have_best = false;
        double min_val = 0;

        for(const auto& pairs : restricts_vector) {
            double val;
            if(restricts_type == 1) {
                val = our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            }
            else {
                val = GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            }
            if(have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }
    }

    if( use_in_exclusion_optimize) {
        int* I = new int[size];
        for(int i = 0; i < size; ++i) I[i] = i;
        in_exclusion_optimize_num = get_vec_optimize_num(I);
        delete[] I;
        if (in_exclusion_optimize_num <= 1) {
            printf("Can not use in_exclusion_optimize with this schedule\n");
            in_exclusion_optimize_num = 0;
        }
        else {
            printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            init_in_exclusion_optimize();
        }
    }
    else {
        in_exclusion_optimize_num = 0;
    }

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }

    build_loop_invariant();
    if( restricts_type != 0) add_restrict(best_pairs);
    
    set_in_exclusion_optimize_redundancy(complete, max_degree);
}

Schedule::Schedule(const int* _adj_mat, int _size, Graph& complete, uint64_t _max_degree) {
    size = _size;
    adj_mat = new int[size * size];
    max_degree = _max_degree;

    memcpy(adj_mat, _adj_mat, size * size * sizeof(int));

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i)
    {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            printf("invalid schedule!\n");
            assert(0);
        }
    }

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy(complete, max_degree);
}

Schedule::~Schedule() {
    adj_mat = nullptr;
    father_prefix_id = nullptr;
    last = nullptr;
    next = nullptr;
    loop_set_prefix_id = nullptr;
    prefix = nullptr;
    restrict_last = nullptr;
    restrict_next = nullptr;
    restrict_index = nullptr;

    // delete[] adj_mat;
    // delete[] father_prefix_id;
    // delete[] last;
    // delete[] next;
    // delete[] loop_set_prefix_id;
    // delete[] prefix;
    // delete[] restrict_last;
    // delete[] restrict_next;
    // delete[] restrict_index;
}

int Schedule::get_max_degree() const {
    int mx = 0;
    for(int i = 0; i < size; ++i) {
        int cnt = 0;
        for(int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i,j,size)];
        if(cnt > mx) mx = cnt;
    }
    return mx;
}

void Schedule::build_loop_invariant() {
    int* tmp_data = new int[size];
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i)
    {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
                tmp_data[data_size++] = j;
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    assert(total_prefix_num <= size * (size - 1) / 2);
    delete[] tmp_data;
}

int Schedule::find_father_prefix(int data_size, const int* data) {
    if (data_size == 0)
        return -1;
    int num = data[data_size - 1];
    for (int prefix_id = last[num]; prefix_id != -1; prefix_id = next[prefix_id])
        if (prefix[prefix_id].equal(data_size, data))
            return prefix_id;
    
    // not found, create new prefix and find its father prefix id recursively
    int father = find_father_prefix(data_size - 1, data);
    father_prefix_id[total_prefix_num] = father;
    next[total_prefix_num] = last[num];
    last[num] = total_prefix_num;
    prefix[total_prefix_num].init(data_size, data);
    ++total_prefix_num;
    return total_prefix_num - 1;
}

void Schedule::add_restrict(int* restrict_first, int* restrict_second, int restrict_size) {
    rest_first = new int[size-1];
    rest_second = new int[size-1];
    memcpy(rest_first, restrict_first, sizeof(int)*(size-1));
    memcpy(rest_second, restrict_second, sizeof(int)*(size-1));
    for(unsigned int i = 0; i < restrict_pair.size(); ) {
        bool tag = true;
        for(unsigned int j = 0; j < restrict_pair.size(); ++j) {
            if(i != j && restrict_pair[j].first == restrict_pair[i].first) 
                for(unsigned int k = 0; k < restrict_pair.size(); ++k)
                    if( i != k && j != k && restrict_pair[k].second == restrict_pair[i].second && restrict_pair[j].second == restrict_pair[k].first ) {
                        tag = false;
                        break;
                    }
            if(tag == false) break;
        }
        if(tag == false) {
            restrict_pair.erase(restrict_pair.begin() + i);
        }
        else ++i;
    }

    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    for (const auto& p : restrict_pair) {
        // p.first must be greater than p.second
        restrict_index[total_restrict_num] = p.first;
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] = total_restrict_num;
        ++total_restrict_num;
    }
}

int Schedule::get_multiplicity() const {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int** isomorphism_vec = new int*[pow];
    for(int i=0; i<pow; i++) isomorphism_vec[i] = new int[size];

    int iso_size = get_isomorphism_vec(isomorphism_vec);

    for(int i=0;i<pow;i++) delete[] isomorphism_vec[i];
    delete[] isomorphism_vec;
    return iso_size;
}

// Schedule::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
int Schedule::aggressive_optimize_get_all_pairs(int** restricts_first, int** restricts_second, int* restricts_size) {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int** isomorphism_vec = new int*[pow];
    for(int i=0; i<pow; i++) isomorphism_vec[i] = new int[size];
    int iso_size = get_isomorphism_vec(isomorphism_vec);

    std::vector<std::vector<std::vector<int>>> permutation_groups;
    permutation_groups.clear();
    for (int i=0; i<iso_size; i++) {
        int* v = new int[size];
        memcpy(v, isomorphism_vec[i], sizeof(int)*size);
        permutation_groups.push_back(calc_permutation_group(v, size));
    }

    // std::vector< std::pair<int,int> > ordered_pairs;
    // ordered_pairs.clear();
    int* ordered_first = new int[permutation_groups.size()];
    int* ordered_second = new int[permutation_groups.size()];
    int ordered_size = 0;

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size(); ) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int>& v : permutation_groups[i]) {
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
            }
            else if (v.size() != 1) {
                two_element_number = -1;
                break;
            }
        }
        if (two_element_number == 1) {
            permutation_groups.erase(permutation_groups.begin() + i);
            for (int j=i+1; j<iso_size; j++) {
                memcpy(isomorphism_vec[j-1], isomorphism_vec[j], sizeof(int)*size);
            }
            iso_size--;
            ordered_first[ordered_size] = found_pair.first;
            ordered_second[ordered_size] = found_pair.second;
            ordered_size++;
            assert(found_pair.first < found_pair.second);
        }
        else
            ++i;
    }

    Pattern base_dag(size);
    for (int i=0; i<ordered_size; i++) {
        base_dag.add_ordered_edge(ordered_first[i], ordered_second[i]);
    }
    
    int rest_size = 0;
    aggressive_optimize_dfs(base_dag, isomorphism_vec, iso_size, permutation_groups, ordered_first, ordered_second, ordered_size, restricts_first, restricts_second, restricts_size, rest_size);

    delete[] ordered_first;
    delete[] ordered_second;
    for(int i=0;i<pow;i++) delete[] isomorphism_vec[i];
    delete[] isomorphism_vec;
    return rest_size;
}

void Schedule::aggressive_optimize_dfs(Pattern base_dag, int** isomorphism_vec, int iso_size, std::vector<std::vector<std::vector<int>>> permutation_groups, int* ordered_first, int* ordered_second, int ordered_size, int** restricts_first, int** restricts_second, int* restricts_size, int& rest_size) {
    for (int i = 0; i < iso_size; ) {
        Pattern test_dag(base_dag);
        int* iso = new int[size];
        memcpy(iso, isomorphism_vec[i], sizeof(int)*size);
        for (int j=0; j<ordered_size; j++) {
            test_dag.add_ordered_edge(iso[ordered_first[j]], iso[ordered_second[j]]);
        }
        delete[] iso;
        if (test_dag.is_dag() == false) {
            permutation_groups.erase(permutation_groups.begin() + i);
            for (int j=i+1; j<iso_size; j++) {
                memcpy(isomorphism_vec[j-1], isomorphism_vec[j], sizeof(int)*size);
            }
            iso_size--;
        }
        else
            ++i;
    }

    if(iso_size == 1) {
        memcpy(restricts_first[rest_size], ordered_first, sizeof(int)*ordered_size);
        memcpy(restricts_second[rest_size], ordered_second, sizeof(int)*ordered_size);
        restricts_size[rest_size] = ordered_size;
        rest_size++;
    }

    std::pair<int, int> found_pair;
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
                std::vector<std::vector< std::vector<int>>> next_permutation_groups = permutation_groups;

                int** next_isomorphism_vec = new int*[iso_size];
                for(int j=0; j<iso_size; j++) next_isomorphism_vec[j] = new int[size];
                memcpy(next_isomorphism_vec, isomorphism_vec, sizeof(int)*size*iso_size);
                int iso_size_t = iso_size;

                int* next_ordered_first = new int[size];
                int* next_ordered_second = new int[size];
                memcpy(next_ordered_first, ordered_first, sizeof(int)*(ordered_size+1));
                memcpy(next_ordered_second, ordered_second, sizeof(int)*(ordered_size+1));
                int next_ordered_size = ordered_size;

                Pattern next_base_dag = base_dag;
                
                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                for (int j=i+1; j<iso_size; j++) {
                    memcpy(next_isomorphism_vec[j-1], next_isomorphism_vec[j], sizeof(int)*size);
                }
                iso_size_t--;
                assert(found_pair.first < found_pair.second);
                next_ordered_first[ordered_size] = found_pair.first;
                next_ordered_second[ordered_size] = found_pair.second;
                next_ordered_size++;
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);

                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, iso_size_t, next_permutation_groups, next_ordered_first, next_ordered_second, ordered_size, restricts_first, restricts_second, restricts_size, rest_size);
                for(int j=0; j<iso_size; j++) delete[] next_isomorphism_vec[j];
                delete[] next_isomorphism_vec;
            }
        if( two_element_number >= 1) {
            break;
        }
        else {
           ++i;
        }
    }
}

void Schedule::GraphZero_get_automorphisms(std::vector<std::vector<int>> &Aut) const {
    int p[size];
    Aut.clear();
    for (int i = 0; i < size; ++i) p[i] = i;
    do {
        bool tag = true;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (adj_mat[INDEX(i, j, size)] != adj_mat[INDEX(p[i], p[j], size)]) {
                    tag = false;
                    break;
                }
            }
            if (!tag) break;
        }
        if (tag) {
            std::vector<int> tmp;
            tmp.clear();
            for (int i = 0; i < size; ++i) tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while (std::next_permutation(p, p + size));
}

int Schedule::get_isomorphism_vec(int** isomorphism_vec) const {
    int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    int** vec = new int*[pow];
    for(int i=0; i<pow; i++) vec[i] = new int[size];
    get_full_permutation(vec);

    int iso_size = 0;
    for(int i=0; i<pow; i++) {
        int* v = new int[size];
        memcpy(v, vec[i], sizeof(int)*size);
        bool flag = true;
        for (int i = 0; i < size; ++i)
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != 0)
                    if (adj_mat[INDEX(v[i], v[j], size)] == 0) // not isomorphism
                    {
                        flag = false;
                        break;
                    }
        if (flag == true) {
            memcpy(isomorphism_vec[iso_size], v, sizeof(int)*size);
            iso_size++;
        }
    }
    return iso_size;
}

void Schedule::get_full_permutation(int** vec) const {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int* tmp = new int[size];
    for (int i=0; i<size; i++) tmp[i] = i;
    memcpy(vec[0], tmp, sizeof(int)*size);
    for (int i=1; i<pow; i++) {
        std::next_permutation(tmp, tmp+size);
        memcpy(vec[i], tmp, sizeof(int)*size);
    }
    delete[] tmp;
}

// ??????????????? like: (A)(B,C)(D)
std::vector<std::vector<int>> Schedule::calc_permutation_group(int* vec, int size) {
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<std::vector<int>> res;
    res.clear();
    for (int i = 0; i < size; ++i) {
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

void Schedule::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;
    
    assert( in_exclusion_optimize_num > 1);

    int* id;
    id = new int[ optimize_num ];

    int* in_exclusion_val;
    in_exclusion_val = new int[ optimize_num * 2];

    for(int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[ 2 * n - 2 ] = 0;
        in_exclusion_val[ 2 * n - 1 ] = 0;

        if( n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int,int> edge[m];
        int e_cnt = 0;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i,j);

        for(int s = 0; s < (1<<m); ++s) {
            dsu.init();
            int bit_cnt = 0;
            for(int i = 0; i < m; ++i) 
                if( s & (1<<i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
            if( dsu.get_set_size() == 1) {
                if( bit_cnt & 1) ++in_exclusion_val[2 * n -1];
                else ++in_exclusion_val[ 2 * n - 2];
            }
        }
    }        

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);

    delete[] id;
    delete[] in_exclusion_val;
}

void Schedule::get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val) {
    if( depth == in_exclusion_optimize_num) {
        int* size = new int[id_cnt];
        for(int i = 0; i < id_cnt; ++i)
            size[i] = 0;
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            size[ id[i] ] ++;
        int val[2];
        val[0] = in_exclusion_val[ size[0] * 2 - 2 ];
        val[1] = in_exclusion_val[ size[0] * 2 - 1 ];
        for(int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];

            val[0] = tmp0 * in_exclusion_val[ size[i] * 2 - 2] + tmp1 * in_exclusion_val[ size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[ size[i] * 2 - 1] + tmp1 * in_exclusion_val[ size[i] * 2 - 2];
        }

        std::vector< std::vector<int> > group;
        group.clear();
        for(int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for(int j = 0; j < in_exclusion_optimize_num; ++j)
                if( id[j] == i) cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back( val[0] - val[1] );

        delete[] size;
        return;
    }
    
    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);
    
    for(int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

int Schedule::restricts_generate(const int* cur_adj_mat, int** restricts_first, int** restricts_second, int* restricts_size, Graph& complete, uint64_t max_degree) {
    Schedule schedule(cur_adj_mat, get_size(), complete, max_degree);
    int rest_size = schedule.aggressive_optimize_get_all_pairs(restricts_first, restricts_second, restricts_size);
    int m = schedule.get_multiplicity();
    long long ans = pattern_matching(complete, schedule, max_degree) / m;
    for (int i = 0; i < rest_size;) {
        Schedule cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size(), complete, max_degree);
        cur_schedule.add_restrict(restricts_first[i], restricts_second[i], restricts_size[i]);
        long long cur_ans = pattern_matching(complete, cur_schedule, max_degree);
        if (cur_ans != ans) {
            for (int j=i+1; j<restricts_size[i]; j++) {
                memcpy(restricts_first[j-1], restricts_first[j], sizeof(int)*(size-1));
                memcpy(restricts_second[j-1], restricts_second[j], sizeof(int)*(size-1));
                restricts_size[j-1] = restricts_size[j];
            }
            rest_size--;
        }
        else {
            ++i;
        }
    }
    return rest_size;
}

int Schedule::get_vec_optimize_num(int* vec) {
    bool is_valid = true;
    for(int i = 1; i < size; ++i) {
        bool have_edge = false;
        for(int j = 0; j < i; ++j)
            if (adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
            }
        if( have_edge == false) {
            is_valid = false;
            break;
        }
    }
    if (!is_valid) return -1;

    for(int k = 2; k <= size; ++k) {
        bool flag = true;
        for(int i = size - k + 1; i < size; ++i) {
            if(adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        }
        if(flag == false) return k - 1;
    }
    assert(0);
    return -1;
}

double Schedule::our_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt; 
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int>> restricts(restrict_size);
    for (int i=0; i<restrict_size; i++) restricts.push_back({restrict_first[i], restrict_second[i]});
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();
    
    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for(int j = 0; j < (int)invariant_size[i].size(); ++j)
            if(invariant_size[i][j] > 1) 
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        if( i ) {
            val *= p_size[1] * pp_size[ cnt_forward - 1 ];
        }
        else {
            val *= p_size[0];
        }
    }
    delete[] cur_adj_mat;

    return val;
}

double Schedule::GraphZero_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt) {
    int max_degree = get_max_degree();
    
    double p_size[max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int>> restricts(restrict_size);
    for (int i=0; i<restrict_size; i++) restricts.push_back({restrict_first[i], restrict_second[i]});
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();
    
    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for(int j = 0; j < (int)invariant_size[i].size(); ++j)
            if(invariant_size[i][j] > 1) 
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= p_size[cnt_forward];

    }
    
    delete[] cur_adj_mat;

    return val;
}

double Schedule::Naive_estimate_schedule_restrict(int* order, int* restrict_first, int* restrict_second, int restrict_size, int v_cnt, unsigned int e_cnt) {
    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];
    
    std::vector<std::pair<int, int>> restricts(restrict_size);
    for (int i=0; i<restrict_size; i++) restricts.push_back({restrict_first[i], restrict_second[i]});
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    int tmp[size];
    
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= v_cnt;
        for(int j = 0; j < i - cnt_forward; ++j)
            val *= (1-p);
        for(int j = 0; j < cnt_forward; ++j)
            val *= p;
    }

    delete[] cur_adj_mat;

    return val;
}

// the i's vertex must connect to at lease one of i-1 vertices at front
int Schedule::remove_invalid_permutation(int** candidate_permutations, int n) {
    int perm_size = n;
    for (int i = 0; i < n; ) {
        const auto& vec = candidate_permutations[i];
        bool tag = true;
        for(int x = 1; x < size; ++x) {
            bool have_edge = false;
            for(int y = 0; y < x; ++y)
                if(adj_mat[INDEX(vec[x],vec[y],size)]) {
                    have_edge = true;
                    break;
                }
            if(!have_edge) {
                tag = false;
                break;
            }
        }
        if(tag) {
            ++i;
        }
        else {
            for (int j=i+1; j<n; j++) memcpy(candidate_permutations[j-1], candidate_permutations[j], sizeof(int)*size);
            perm_size--;
        }
    }
    return perm_size;
}

void Schedule::set_in_exclusion_optimize_redundancy(Graph& complete, uint64_t max_degree) {
    int tmp = get_in_exclusion_optimize_num();
    if (tmp <= 1) {
        in_exclusion_optimize_redundancy = 1;
    }
    else {
        in_exclusion_optimize_redundancy = 1;
        long long ans = pattern_matching(complete, *this, max_degree);
        set_in_exclusion_optimize_num(0);
        long long true_ans = pattern_matching(complete, *this, max_degree);
        set_in_exclusion_optimize_num(tmp);
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}
