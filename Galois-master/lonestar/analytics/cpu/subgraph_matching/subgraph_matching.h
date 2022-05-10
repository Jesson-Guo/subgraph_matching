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
#include <mutex>

std::mutex vertex_mutex;

void quickSort(int* first, int* second, int l, int r) {
    if (l < r) {      
        int i = l, j = r, x = first[l], y = second[l];
        while (i < j) {
            while(i < j && (first[j] > x || (first[j] == x && second[j] >= y))) j--;
            if(i < j) {
                first[i] = first[j];
                second[i] = second[j];
                i++;
            }
            while(i < j && first[i] < x) i++;
            if(i < j) {
                first[j] = first[i];
                second[j] = second[i];
                j--;
            }
        }
        first[i] = x;
        second[i] = y;
        quickSort(first, second, l, i - 1);
        quickSort(first, second, i + 1, r);
    }
}

void build_vertex_set_with_block(const Schedule& schedule, VertexSet* vertex_set, int* input_data, int input_size, int prefix_id) {
    std::unique_lock<std::mutex> lk(vertex_mutex);
    vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, input_data, input_size, prefix_id);
}

void push_back_mutex(VertexSet& set, int v) {
    std::unique_lock<std::mutex> lk(vertex_mutex);
    set.push_back(v);
}

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
                // 这里可以输出具体的匹配结果
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
            build_vertex_set_with_block(schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id);
            // vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id, vertex);
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
    // 给它们加锁
    VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
    galois::GAccumulator<long long> local_ans;
    // TODO : try different chunksize
    galois::do_all(
        galois::iterate(graph),
        [&](uint64_t vertex) {
            VertexSet subtraction_set;
            VertexSet tmp_set;
            subtraction_set.init();
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                build_vertex_set_with_block(
                    schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id);
                // vertex_set[prefix_id].build_vertex_set(
                //     schedule, vertex_set, &graph.getData(vertex).adj_node[0], graph.getData(vertex).degree, prefix_id);
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
    if (performance_modeling_type != 0 && tri_cnt == -1) {
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

    int* best_pairs_first = new int[size-1];
    int* best_pairs_second = new int[size-1];
    int best_pairs_size = 0;
    //Initialize adj_mat
    //If we use performance_modeling, we may change the order of vertex,
    //the best order produced by performance_modeling(...) is saved in best_order[]
    //Finally, we use best_order[] to relocate adj_mat
    if (performance_modeling_type != 0) { 
        int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;

        int** candidate_permutations = new int*[pow];
        for (int i=0; i<pow; i++) candidate_permutations[i] = new int[size];

        get_full_permutation(candidate_permutations);

        int perm_size = remove_invalid_permutation(candidate_permutations, pow);
        int** candidate_permutations_t = new int*[pow];
        for (int i=0; i<pow; i++) candidate_permutations_t[i] = new int[size];

        if (performance_modeling_type == 1) {
            //reduce candidates
            int max_val = 0;
            for (int i=0; i<perm_size; i++) {
                max_val = std::max(max_val, get_vec_optimize_num(candidate_permutations[i]));
            }
            int j = 0;
            for (int i=0; i<perm_size; i++) {
                if (get_vec_optimize_num(candidate_permutations[i]) == max_val) {
                    memcpy(candidate_permutations_t[j], candidate_permutations[i], sizeof(int)*size);
                    j++;
                }
            }
            perm_size = j;
        }
        for (int i=0;i<pow;i++) delete[] candidate_permutations[i];
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

        for (int o=0; o<perm_size; o++) {
            int* rank = new int[size];
            // int rank[size];
            for (int i = 0; i < size; ++i) rank[candidate_permutations_t[o][i]] = i;

            int* cur_adj_mat = new int[size*size];
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
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

            if (restricts_type == 1) {
                rest_size = restricts_generate(cur_adj_mat, restricts_first, restricts_second, restricts_size, complete, max_degree);
            }

            if (rest_size == 0) {
                double val = 0;
                if (performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(candidate_permutations_t[o], nullptr, nullptr, 0, v_cnt, e_cnt, tri_cnt);
                }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = candidate_permutations_t[o][i];
                    best_pairs_size = 0;
                }
            }

            // for (const auto& pairs : restricts_vector) {
            for (int i=0; i<rest_size; i++) {
                double val = 0;
                if (performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(candidate_permutations_t[o], restricts_first[i], restricts_second[i], restricts_size[i], v_cnt, e_cnt, tri_cnt);
                }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = candidate_permutations_t[o][i];
                    best_pairs_size = restricts_size[i];
                    memcpy(best_pairs_first, restricts_first[i], sizeof(int)*restricts_size[i]);
                    memcpy(best_pairs_second, restricts_second[i], sizeof(int)*restricts_size[i]);
                }
            }
            // delete[] rank;
            // delete[] cur_adj_mat;
            // delete[] restricts_size;
            // for (int i=0;i<pow;i++) delete[] restricts_first[i];
            // for (int i=0;i<pow;i++) delete[] restricts_second[i];
            // delete[] restricts_first;
            // delete[] restricts_second;
        }
        for (int i=0;i<pow;i++) delete[] candidate_permutations_t[i];
        delete[] candidate_permutations_t;

        int rank[size];
        for (int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const int* pattern_adj_mat = pattern.get_adj_mat_ptr();
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)]; 
        delete[] best_order;
    }
    if (use_in_exclusion_optimize) {
        int* I = new int[size];
        for (int i = 0; i < size; ++i) I[i] = i;
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
        if (valid == false) {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }

    build_loop_invariant();
    if (restricts_type != 0) add_restrict(best_pairs_first, best_pairs_second, best_pairs_size);
    delete[] best_pairs_first;
    delete[] best_pairs_second;
    
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
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j) {
            if (adj_mat[INDEX(i, j, size)]) {
                valid = true;
                break;
            }
        }
        if (valid == false) {
            printf("invalid schedule!\n");
            assert(0);
        }
    }

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy(complete, max_degree);
}

Schedule::~Schedule() {
    // delete[] adj_mat;
    // delete[] father_prefix_id;
    // delete[] last;
    // delete[] next;
    // delete[] loop_set_prefix_id;
    // delete[] prefix;
    // delete[] restrict_last;
    // delete[] restrict_next;
    // delete[] restrict_index;

    // delete[] rest_first;
    // delete[] rest_second;
}

int Schedule::get_max_degree() const {
    int mx = 0;
    for (int i = 0; i < size; ++i) {
        int cnt = 0;
        for (int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i,j,size)];
        if (cnt > mx) mx = cnt;
    }
    return mx;
}

void Schedule::build_loop_invariant() {
    int tmp_data[size];
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
    rest_size = restrict_size;
    memcpy(rest_first, restrict_first, sizeof(int)*(size-1));
    memcpy(rest_second, restrict_second, sizeof(int)*(size-1));
    for (int i = 0; i < rest_size; ) {
        bool tag = true;
        for (int j = 0; j < rest_size; ++j) {
            if (i != j && rest_first[j] == rest_first[i]) { 
                for (int k = 0; k < rest_size; ++k) {
                    if (i != k && j != k && rest_second[k] == rest_second[i] && rest_second[j] == rest_first[k]) {
                        tag = false;
                        break;
                    }
                }
            }
            if (tag == false) break;
        }
        if (tag == false) {
            for (int j=i+1; j<rest_size; j++) {
                rest_first[i-1] = rest_first[i];
                rest_second[i-1] = rest_second[i];
            }
            rest_size--;
        }
        else ++i;
    }

    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    for (int i=0; i<rest_size; i++) {
        restrict_index[total_restrict_num] = rest_first[i];
        restrict_next[total_restrict_num] = restrict_last[rest_second[i]];
        restrict_last[rest_second[i]] = total_restrict_num;
        ++total_restrict_num;
    }
}

int Schedule::get_multiplicity() const {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int** isomorphism_vec = new int*[pow];
    for (int i=0; i<pow; i++) isomorphism_vec[i] = new int[size];

    int iso_size = get_isomorphism_vec(isomorphism_vec);

    // for (int i=0;i<pow;i++) delete[] isomorphism_vec[i];
    // delete[] isomorphism_vec;
    return iso_size;
}

// Schedule::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
int Schedule::aggressive_optimize_get_all_pairs(int** restricts_first, int** restricts_second, int* restricts_size) {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int** isomorphism_vec = new int*[pow];
    for (int i=0; i<pow; i++) isomorphism_vec[i] = new int[size];
    int iso_size = get_isomorphism_vec(isomorphism_vec);

    int*** permutation_groups = new int**[pow];
    for (int i = 0;i < pow;i++) {
        permutation_groups[i] = new int*[size];
        for (int j = 0;j < size;j++) {
            permutation_groups[i][j] = new int[size];
        }
    }
    int** permutation_groups_size = new int*[pow];
    for (int j = 0;j < pow;j++) {
        permutation_groups_size[j] = new int[size];
    }
    int groups_cnt[pow];

    for (int i=0; i<iso_size; i++) {
        groups_cnt[i] = calc_permutation_group(isomorphism_vec[i], size, permutation_groups[i], permutation_groups_size[i]);
    }

    int ordered_first[iso_size];
    int ordered_second[iso_size];
    int ordered_size = 0;

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (int i = 0; i < iso_size; ) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (int j = 0; j < groups_cnt[i]; j++) {
            if (permutation_groups_size[i][j] == 2) {
                ++two_element_number;
                found_pair = std::pair<int ,int>(permutation_groups[i][j][0], permutation_groups[i][j][1]);
            }
            else if (permutation_groups_size[i][j] != 1) {
                two_element_number = -1;
                break;
            }
        }
        if (two_element_number == 1) {
            for (int j=i+1; j<iso_size; j++) {
                memcpy(isomorphism_vec[j-1], isomorphism_vec[j], sizeof(int)*size);
                memcpy(permutation_groups[j-1], permutation_groups[j], sizeof(int)*size*size);
                memcpy(permutation_groups_size[j-1], permutation_groups_size[j], sizeof(int)*size);
                groups_cnt[j-1] = groups_cnt[j];
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
    aggressive_optimize_dfs(base_dag, isomorphism_vec, iso_size, permutation_groups, permutation_groups_size, groups_cnt, ordered_first, ordered_second, ordered_size, restricts_first, restricts_second, restricts_size, rest_size);

    // for (int i = 0;i < pow;i++) {
    //     for (int j = 0;j < size;j++) delete[] permutation_groups[i][j];
    //     delete[] permutation_groups[i];
    // }
    // delete[] permutation_groups;
    // for (int j = 0;j < pow;j++) delete[] permutation_groups_size[j];
    // delete[] permutation_groups_size;
    // delete[] groups_cnt;

    // for (int i=0;i<pow;i++) delete[] isomorphism_vec[i];
    // delete[] isomorphism_vec;
    return rest_size;
}

void Schedule::aggressive_optimize_dfs(Pattern base_dag, int** isomorphism_vec, int iso_size, int*** permutation_groups, int** permutation_groups_size, int* groups_cnt, int* ordered_first, int* ordered_second, int ordered_size, int** restricts_first, int** restricts_second, int* restricts_size, int& rest_size) {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    for (int i = 0; i < iso_size; ) {
        Pattern test_dag(base_dag);
        for (int j=0; j<ordered_size; j++) {
            test_dag.add_ordered_edge(isomorphism_vec[i][ordered_first[j]], isomorphism_vec[i][ordered_second[j]]);
        }
        if (test_dag.is_dag() == false) {
            for (int j=i+1; j<iso_size; j++) {
                memcpy(isomorphism_vec[j-1], isomorphism_vec[j], sizeof(int)*size);
                memcpy(permutation_groups[j-1], permutation_groups[j], sizeof(int)*size*size);
                memcpy(permutation_groups_size[j-1], permutation_groups_size[j], sizeof(int)*size);
                groups_cnt[j-1] = groups_cnt[j];
            }
            iso_size--;
        }
        else
            ++i;
    }

    if (iso_size == 1) {
        memcpy(restricts_first[rest_size], ordered_first, sizeof(int)*ordered_size);
        memcpy(restricts_second[rest_size], ordered_second, sizeof(int)*ordered_size);
        restricts_size[rest_size] = ordered_size;
        rest_size++;
    }

    std::pair<int, int> found_pair;
    for (int i = 0; i < iso_size; ) {
        int two_element_number = 0;
        for (int j = 0; j < groups_cnt[i]; j++) {
            if (permutation_groups_size[i][j] == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(permutation_groups[i][j][0], permutation_groups[i][j][1]);

                int*** next_permutation_groups = new int**[pow];
                for (int i = 0;i < pow;i++) {
                    next_permutation_groups[i] = new int*[size];
                    for (int j = 0;j < size;j++) {
                        next_permutation_groups[i][j] = new int[size];
                    }
                }
                int** next_permutation_groups_size = new int*[pow];
                for (int j = 0;j < pow;j++) {
                    next_permutation_groups_size[j] = new int[size];
                }
                int next_groups_cnt[pow];
                memcpy(next_permutation_groups, permutation_groups, sizeof(int)*size*size*pow);
                memcpy(next_permutation_groups_size, permutation_groups_size, sizeof(int)*size*pow);
                memcpy(next_groups_cnt, groups_cnt, sizeof(int)*pow);

                int** next_isomorphism_vec = new int*[pow];
                for (int j=0; j<pow; j++) next_isomorphism_vec[j] = new int[size];
                memcpy(next_isomorphism_vec, isomorphism_vec, sizeof(int)*size*pow);
                int iso_size_t = iso_size;

                int next_ordered_first[size];
                int next_ordered_second[size];
                memcpy(next_ordered_first, ordered_first, sizeof(int)*(ordered_size+1));
                memcpy(next_ordered_second, ordered_second, sizeof(int)*(ordered_size+1));
                int next_ordered_size = ordered_size;

                Pattern next_base_dag = base_dag;

                for (int j=i+1; j<iso_size; j++) {
                    memcpy(next_isomorphism_vec[j-1], next_isomorphism_vec[j], sizeof(int)*size);
                    memcpy(next_permutation_groups[j-1], next_permutation_groups[j], sizeof(int)*size*size);
                    memcpy(next_permutation_groups_size[j-1], next_permutation_groups_size[j], sizeof(int)*size);
                    groups_cnt[j-1] = groups_cnt[j];
                }
                iso_size_t--;
                assert(found_pair.first < found_pair.second);
                next_ordered_first[ordered_size] = found_pair.first;
                next_ordered_second[ordered_size] = found_pair.second;
                next_ordered_size++;
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);

                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, iso_size_t, next_permutation_groups, next_permutation_groups_size, next_groups_cnt, next_ordered_first, next_ordered_second, ordered_size, restricts_first, restricts_second, restricts_size, rest_size);

                // for (int i = 0;i < pow;i++) {
                //     for (int j = 0;j < size;j++) delete[] next_permutation_groups[i][j];
                //     delete[] next_permutation_groups[i];
                // }
                // delete[] next_permutation_groups;
                // for (int j = 0;j < pow;j++) delete[] next_permutation_groups_size[j];
                // delete[] next_permutation_groups_size;

                // for (int j=0; j<pow; j++) delete[] next_isomorphism_vec[j];
                // delete[] next_isomorphism_vec;
            }
        }
        if (two_element_number >= 1) {
            break;
        }
        else {
           ++i;
        }
    }
}

int Schedule::get_isomorphism_vec(int** isomorphism_vec) const {
    int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    int** vec = new int*[pow];
    for (int i=0; i<pow; i++) vec[i] = new int[size];
    get_full_permutation(vec);

    int iso_size = 0;
    for (int i=0; i<pow; i++) {
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
        // delete[] v;
    }

    // for (int j=0; j<iso_size; j++) delete[] vec[j];
    // delete[] vec;
    return iso_size;
}

void Schedule::get_full_permutation(int** vec) const {
    int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;
    int tmp[size];
    for (int i=0; i<size; i++) tmp[i] = i;
    memcpy(vec[0], tmp, sizeof(int)*size);
    for (int i=1; i<pow; i++) {
        std::next_permutation(tmp, tmp+size);
        memcpy(vec[i], tmp, sizeof(int)*size);
    }
}

// 计算置换群 like: (A)(B,C)(D)
int Schedule::calc_permutation_group(int* vec, int size, int** perm_group, int* perm_group_size) {
    bool use[size];
    for (int i = 0; i < size; ++i) use[i] = false;
    int j=0;
    for (int i = 0; i < size; ++i) {
        if (use[i] == false) {
            perm_group[j][perm_group_size[j]] = i;
            perm_group_size[j]++;
            use[i] = true;
            int x = vec[i];
            while (use[x] == false) {
                use[x] = true;
                perm_group[j][perm_group_size[j]] = x;
                perm_group_size[j]++;
                x = vec[x];
            }
            j++;
        }
    }
    return j;
}

void Schedule::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;
    assert( in_exclusion_optimize_num > 1);
    int id[optimize_num];
    int in_exclusion_val[optimize_num * 2];

    for (int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[2 * n - 2 ] = 0;
        in_exclusion_val[2 * n - 1 ] = 0;

        if (n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int,int> edge[m];
        int e_cnt = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i,j);

        for (int s = 0; s < (1<<m); ++s) {
            dsu.init();
            int bit_cnt = 0;
            for (int i = 0; i < m; ++i) 
                if (s & (1<<i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
            if (dsu.get_set_size() == 1) {
                if (bit_cnt & 1) ++in_exclusion_val[2 * n -1];
                else ++in_exclusion_val[ 2 * n - 2];
            }
        }
    }        

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);
}

void Schedule::get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val) {
    if (depth == in_exclusion_optimize_num) {
        int opt_size[id_cnt];
        for (int i = 0; i < id_cnt; ++i) opt_size[i] = 0;
        for (int i = 0; i < in_exclusion_optimize_num; ++i) opt_size[id[i]] ++;
        int val[2];
        val[0] = in_exclusion_val[opt_size[0] * 2 - 2];
        val[1] = in_exclusion_val[opt_size[0] * 2 - 1];
        for (int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];
            val[0] = tmp0 * in_exclusion_val[opt_size[i] * 2 - 2] + tmp1 * in_exclusion_val[opt_size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[opt_size[i] * 2 - 1] + tmp1 * in_exclusion_val[opt_size[i] * 2 - 2];
        }

        std::vector<std::vector<int>> group;
        group.clear();
        for (int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for (int j = 0; j < in_exclusion_optimize_num; ++j)
                if (id[j] == i) cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back( val[0] - val[1] );
        return;
    }
    
    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);
    
    for (int i = 0; i < id_cnt; ++i) {
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
    for (int i = 1; i < size; ++i) {
        bool have_edge = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
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
            if (adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        }
        if (flag == false) return k - 1;
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
    for (int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;
    
    int cur_adj_mat[size*size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    quickSort(restrict_first, restrict_second, 0, restrict_size-1);

    double sum[restrict_size];
    for (int i = 0; i < restrict_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restrict_size; ++i)
            if (tmp[restrict_first[i]] > tmp[restrict_second[i]]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;
    for (int i = 0; i < restrict_size; ++i)
        sum[i] = sum[i] /total;
    for (int i = restrict_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    int** invariant_size = new int*[size];
    for (int i=0; i<size; i++) invariant_size[i] = new int[size];
    int* invariant_size_cnt = new int[size];
    for (int i=0; i<size; i++) invariant_size_cnt[i] = 0;
    
    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j) {
            if (cur_adj_mat[INDEX(j, i, size)]) {
                invariant_size[j][invariant_size_cnt[j]] = c--;
                invariant_size_cnt[j]++;
            }
        }

        for (int j = 0; j < invariant_size_cnt[i]; ++j)
            if (invariant_size[i][j] > 1) 
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for (int j = 0; j < restrict_size; ++j)
            if (restrict_second[j] == i)
                val *=  sum[j];
        if (i ) {
            val *= p_size[1] * pp_size[ cnt_forward - 1 ];
        }
        else {
            val *= p_size[0];
        }
    }
    return val;
}

// the i's vertex must connect to at lease one of i-1 vertices at front
int Schedule::remove_invalid_permutation(int** candidate_permutations, int n) {
    int perm_size = n;
    for (int i = 0; i < n; ) {
        const auto& vec = candidate_permutations[i];
        bool tag = true;
        for (int x = 1; x < size; ++x) {
            bool have_edge = false;
            for (int y = 0; y < x; ++y)
                if (adj_mat[INDEX(vec[x],vec[y],size)]) {
                    have_edge = true;
                    break;
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
