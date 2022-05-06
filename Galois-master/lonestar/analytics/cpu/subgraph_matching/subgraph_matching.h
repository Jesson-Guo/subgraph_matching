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
static uint32_t max_degree = 0;

//! [withnumaalloc]
using Graph = galois::graphs::LC_CSR_Graph<NodeData, void>::
    with_no_lockable<true>::type ::with_numa_alloc<true>::type;
//! [withnumaalloc]
typedef Graph::GraphNode GNode;

#include "disjoint_set_union.h"
#include "vertex_set.h"
#include "configuration.h"
#include "utils.h"

#include <iostream>
#include <algorithm>

void pattern_matching_aggressive_func(Graph& graph, const Configuration &config, std::vector<VertexSet> vertex_set, VertexSet &subtraction_set, VertexSet &tmp_set, long long &local_ans, int depth) {
    int loop_set_prefix_id = config.get_loop_set_prefix_id(depth); // @@@
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;

    std::vector<int> loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    //Case: in_exclusion_optimize_num > 1
    if (depth == config.get_size() - config.get_in_exclusion_optimize_num()) {
        int in_exclusion_optimize_num = config.get_in_exclusion_optimize_num(); // @@@
        int loop_set_prefix_ids[in_exclusion_optimize_num];
        loop_set_prefix_ids[0] = loop_set_prefix_id;
        for (int i = 1; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = config.get_loop_set_prefix_id(depth + i);
        for (int optimize_rank = 0; optimize_rank < (int)config.in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector<std::vector<int>> &cur_graph = config.in_exclusion_optimize_group[optimize_rank];
            long long val = config.in_exclusion_optimize_val[optimize_rank];
            for (int cur_graph_rank = 0; cur_graph_rank < (int)cur_graph.size(); ++cur_graph_rank) {
                //if size == 1 , we will not call intersection(...)
                //so we will not allocate memory for data
                //otherwise, we need to copy the data to do intersection(...)
                if (cur_graph[cur_graph_rank].size() == 1) {
                    int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    val = val * unorderd_subtraction_size(vertex_set[id], subtraction_set);
                }
                else {
                    int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
                    tmp_set.init(max_degree);
                    tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

                    for (int i = 2; i < (int)cur_graph[cur_graph_rank].size(); ++i) {
                        int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
                        VertexSet tmp_t;
                        tmp_t.init(tmp_set.get_size(), tmp_set.get_data_ptr());
                        tmp_set.intersection(vertex_set[id], tmp_t);
                    }
                    val = val * unorderd_subtraction_size(tmp_set, subtraction_set);
                }
                if (val == 0)
                    break;
            }
            local_ans += val;
        }
        return;
    }
    //Case: in_exclusion_optimize_num <= 1
    if (depth == config.get_size() - 1) {
        // TODO : try more kinds of calculation. @@@
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (config.get_total_restrict_num() > 0) {
            int min_vertex = graph.size();
            for (int i = config.get_restrict_last(depth); i != -1; i = config.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(config.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(config.get_restrict_index(i));
            const VertexSet &vset = vertex_set[loop_set_prefix_id];
            std::vector<int> vset_data = vset.get_data_ptr();
            auto lower = std::lower_bound(vset_data.begin(), vset_data.begin() + vset.get_size(), min_vertex);
            int size_after_restrict = std::distance(vset_data.begin(), lower);
            if (size_after_restrict > 0) {
                // 这里可以输出具体的匹配结果
                local_ans += unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
            }
        }
        else
            local_ans += unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    // TODO : min_vertex is also a loop invariant
    int min_vertex = graph.size();
    for (int i = config.get_restrict_last(depth); i != -1; i = config.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(config.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(config.get_restrict_index(i));
    int ii = 0;
    for (int &i = ii; i < loop_size; ++i) {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        bool is_zero = false;
        for (int prefix_id = config.get_last(depth); prefix_id != -1; prefix_id = config.get_next(prefix_id)) {
            // vertex_set[prefix_id].build_vertex_set(config.get_father_prefix_id(prefix_id), vertex_set, graph.getData(vertex).adj_node, graph.getData(vertex).degree, vertex);
            vertex_set[prefix_id].build_vertex_set(config.get_father_prefix_id(prefix_id), vertex_set, graph.getData(vertex).adj_node, graph.getData(vertex).degree);
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if (is_zero)
            continue;
        subtraction_set.push_back(vertex);
        pattern_matching_aggressive_func(graph, config, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);
        subtraction_set.pop_back();
    }
}

long long pattern_matching(Graph& graph, const Configuration &config) {
    long long global_ans = 0;
    std::vector<VertexSet> vertex_set = std::vector<VertexSet>(config.get_total_prefix_num());
    VertexSet subtraction_set;
    VertexSet tmp_set;
    subtraction_set.init();
    long long local_ans = 0;
    // TODO : try different chunksize
    galois::do_all(
        galois::iterate(graph),
        [&](uint64_t vertex) {
            for (int prefix_id = config.get_last(0); prefix_id != -1; prefix_id = config.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(config.get_father_prefix_id(prefix_id), vertex_set, graph.getData(vertex).adj_node, graph.getData(vertex).degree);
            }
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(graph, config, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            subtraction_set.pop_back();
        },
        galois::loopname("pattern matching"),
        galois::chunk_size<CHUNK_SIZE>(),
        galois::steal(),
        galois::no_stats()
    );
    global_ans += local_ans;
    return global_ans / config.get_in_exclusion_optimize_redundancy();
}

void schedule_generator() {

}

void restricts_generator(const std::string cur_adj_mat, int size, std::vector<std::vector<std::pair<int,int>>> &restricts, Graph& complete) {
    Configuration config(cur_adj_mat, size);
    config.init(complete);
    config.aggressive_optimize_get_all_pairs(restricts);
    long long ans = pattern_matching(complete, config);
    ans /= get_multiplicity(size, cur_adj_mat);
    for (int i = 0; i < (int)restricts.size();) {
        Configuration cur_config(config.get_adj_mat_ptr(), config.get_size());
        cur_config.init(complete);
        cur_config.add_restrict(restricts[i]);
        long long cur_ans = pattern_matching(complete, cur_config);
        if (cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        }
        else ++i;
    }
}

// performance_modeling type = 0 : not use modeling
//                      type = 1 : use our modeling
//                      type = 2 : use GraphZero's modeling
//                      type = 3 : use naive modeling
// restricts_type = 0 : not use restricts
//                = 1 : use our restricts
//                = 2 : use GraphZero's restricts
void configuration_generation(Configuration& config, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize, Graph& complete, int v_cnt, uint32_t e_cnt, long long tri_cnt) {
    if (performance_modeling_type != 0 && tri_cnt == -1) {
        printf("Fatal: Can not use performance modeling if not have triangle number of this dataset.\n");
        fflush(stdout);
        assert(0);
    }

    int size = config.get_size();
    std::string adj_mat = config.get_adj_mat_ptr();

    std::vector<std::pair<int,int>> best_pairs;
    best_pairs.clear();

    if (performance_modeling_type != 0) {
        unsigned int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;
        std::vector<std::vector<int>> candidate_permutations;
        candidate_permutations.clear();

        bool use[size];
        for (int i = 0; i < size; ++i) use[i] = false;
        std::vector<int> tmp_vec;
        get_full_permutation(size, candidate_permutations, use, tmp_vec, 0);
        assert(candidate_permutations.size() == pow);

        remove_invalid_permutation(candidate_permutations, size, adj_mat);

        if (performance_modeling_type == 1) {
            //reduce candidates
            int max_val = 0;
            for (const auto &vec : candidate_permutations) {
                max_val = std::max(max_val, get_vec_optimize_num(size, vec, adj_mat));
            }
            std::vector<std::vector<int>> tmp;
            tmp.clear();
            for (const auto &vec : candidate_permutations) {
                if (get_vec_optimize_num(size, vec, adj_mat) == max_val) {
                    tmp.push_back(vec);
                }
            }
            candidate_permutations = tmp;
        }

        std::vector<int> best_order(size);
        double min_val(std::numeric_limits<double>::max() / 4);
        bool have_best = false;

        for (const auto &vec : candidate_permutations) {
            int rank[size];
            for (int i = 0; i < size; ++i) 
                rank[vec[i]] = i;
            std::string cur_adj_mat(size*size, '0');
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
                    cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

            std::vector<std::vector<std::pair<int,int>>> restricts_vector;

            restricts_vector.clear();

            if (restricts_type == 1) {
                restricts_generator(cur_adj_mat, size, restricts_vector, complete);
            }
            else {
                Configuration conf(cur_adj_mat, size);
                conf.init(complete);
                std::vector<std::pair<int, int>> pairs;
                conf.GraphZero_aggressive_optimize(pairs);
                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }

            if (restricts_vector.size() == 0) {
                std::vector<std::pair<int,int>> Empty;
                Empty.clear();
                double val;
                if (performance_modeling_type == 1) {
                    val = config.our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if (performance_modeling_type == 2) {
                        val = config.GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                    else {
                        val = config.Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                }
                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }

            for (const auto& pairs : restricts_vector) {
                double val;
                if (performance_modeling_type == 1) {
                    val = config.our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if (performance_modeling_type == 2) {
                        val = config.GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                    else {
                        val = config.Naive_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                }
                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }
        }

        int rank[size];
        for (int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const std::string pattern_adj_mat = adj_mat;
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)];
        config.set_adj_mat(adj_mat);
    }
    else {
        std::vector<int> I;
        I.clear();
        for (int i = 0; i < size; ++i) I.push_back(i);

        std::vector<std::vector<std::pair<int,int>>> restricts_vector;
        restricts_vector.clear();

        if (restricts_type != 0) {
            if (restricts_type == 1) {
                restricts_generator(adj_mat, size, restricts_vector, complete);
            }
            else {
                std::vector< std::pair<int,int>> pairs;
                config.GraphZero_aggressive_optimize(pairs);
                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
        }

        bool have_best = false;
        double min_val(std::numeric_limits<double>::max() / 4);

        for (const auto& pairs : restricts_vector) {
            double val;
            if (restricts_type == 1) {
                val = config.our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            }
            else {
                val = config.GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            }
            if (have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }
    }

    if (use_in_exclusion_optimize) {
        std::vector<int> I;
        I.clear();
        for (int i = 0; i < size; ++i) I.push_back(i);
        int in_exclusion_optimize_num = get_vec_optimize_num(size, I, adj_mat);
        if (in_exclusion_optimize_num <= 1) {
            printf("Can not use in_exclusion_optimize with this schedule\n");
            config.set_in_exclusion_optimize_num(0);
        }
        else {
            printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            config.init_in_exclusion_optimize();
        }
    }
    else {
        config.set_in_exclusion_optimize_num(0);
    }

    if (restricts_type != 0) config.add_restrict(best_pairs);
}

// Configuration functions
Configuration::Configuration(const std::string _adj_mat, int _size) {
    adj_mat = _adj_mat;
    size = _size;
    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;

    loop_set_prefix_id = std::vector<int>(size);
    prefix = std::vector<Prefix>(max_prefix_num);
    restrict_index = std::vector<int>(max_prefix_num);

    father_prefix_id = std::vector<int>(max_prefix_num, -1);
    last = std::vector<int>(size, -1);
    next = std::vector<int>(max_prefix_num, -1);
    restrict_last = std::vector<int>(size, -1);
    restrict_next = std::vector<int>(max_prefix_num, -1);

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;
}

void Configuration::init(Graph& complete) {
    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j) {
            if (adj_mat[INDEX(i, j, size)] == '1') {
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

    set_in_exclusion_optimize_redundancy(complete);
}

Configuration::~Configuration() {}

void Configuration::build_loop_invariant() {
    std::vector<int> tmp_data = std::vector<int>(size);
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i) {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)] == '1')
                tmp_data[data_size++] = j;
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    assert(total_prefix_num <= size * (size - 1) / 2);
}

int Configuration::find_father_prefix(int data_size, std::vector<int> data) {
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

void Configuration::add_restrict(const std::vector<std::pair<int, int>>& restricts) {
    restrict_pair = restricts;
    for (unsigned int i = 0; i < restrict_pair.size();) {
        bool tag = true;
        for (unsigned int j = 0; j < restrict_pair.size(); ++j) {
            if (i != j && restrict_pair[j].first == restrict_pair[i].first) {
                for (unsigned int k = 0; k < restrict_pair.size(); ++k) {
                    if (i != k && j != k && restrict_pair[k].second == restrict_pair[i].second && restrict_pair[j].second == restrict_pair[k].first) {
                        tag = false;
                        break;
                    }
                }
            }
            if (tag == false) break;
        }
        if (tag == false) {
            restrict_pair.erase(restrict_pair.begin() + i);
        }
        else ++i;
    }

    int max_prefix_num = size * (size - 1) / 2;
    restrict_last = std::vector<int>(size, -1);
    restrict_next = std::vector<int>(max_prefix_num, -1);
    total_restrict_num = 0;
    for (const auto& p : restrict_pair) {
        // p.first must be greater than p.second
        restrict_index[total_restrict_num] = p.first;
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] = total_restrict_num;
        ++total_restrict_num;
    }
}

// Configuration::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
void Configuration::aggressive_optimize_get_all_pairs(std::vector<std::vector<std::pair<int, int>>>& ordered_pairs_vector) {
    std::vector<std::vector<int>> isomorphism_vec = get_isomorphism_vec(size, adj_mat);

    std::vector<std::vector<std::vector<int>>> permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int>& v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs_vector.clear();

    std::vector<std::pair<int,int>> ordered_pairs;
    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size();) {
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
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        }
        else ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int>& pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    aggressive_optimize_dfs(base_dag, isomorphism_vec, permutation_groups, ordered_pairs, ordered_pairs_vector);
}

void Configuration::aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int>> isomorphism_vec, std::vector<std::vector<std::vector<int>>> permutation_groups, std::vector<std::pair<int,int>> ordered_pairs, std::vector<std::vector<std::pair<int,int>>>& ordered_pairs_vector) {
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
                found_pair = std::pair<int ,int>(v[0], v[1]);
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

void Configuration::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;
    assert(in_exclusion_optimize_num> 1);
    std::vector<int> id(optimize_num);
    std::vector<int> in_exclusion_val(optimize_num*2);

    for (int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[2 * n - 2] = 0;
        in_exclusion_val[2 * n - 1] = 0;

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
            for (int i = 0; i < m; ++i) {
                if (s & (1<<i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
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

void Configuration::get_in_exclusion_optimize_group(int depth, std::vector<int>& id, int id_cnt, std::vector<int>& in_exclusion_val) {
    if (depth == in_exclusion_optimize_num) {
        int* exc_size = new int[id_cnt];
        for (int i = 0; i < id_cnt; ++i) exc_size[i] = 0;
        for (int i = 0; i < in_exclusion_optimize_num; ++i) exc_size[id[i]] ++;
        int val[2];
        val[0] = in_exclusion_val[exc_size[0] * 2 - 2];
        val[1] = in_exclusion_val[exc_size[0] * 2 - 1];
        for (int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];
            val[0] = tmp0 * in_exclusion_val[exc_size[i] * 2 - 2] + tmp1 * in_exclusion_val[exc_size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[exc_size[i] * 2 - 1] + tmp1 * in_exclusion_val[exc_size[i] * 2 - 2];
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
        in_exclusion_optimize_val.push_back(val[0] - val[1]);
        delete[] exc_size;
        return;
    }
    
    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);
    
    for (int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

void Configuration::GraphZero_aggressive_optimize(std::vector<std::pair<int, int>>& ordered_pairs) const { 
    std::vector<std::vector<int>> Aut;
    GraphZero_get_automorphisms(Aut);

    std::vector<std::pair<int, int>> L;
    L.clear();

    for (int v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector<std::vector<int>> stabilized_aut;
        stabilized_aut.clear();

        for (int i = 0; i < (int)Aut.size(); ++i) {
            std::vector<int>& x = Aut[i];
            if (x[v] == v) {
                stabilized_aut.push_back(x);
            }
            else {
                int x1 = v, x2 = x[v];
                if (x1 > x2) {
                    int tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<int, int> cur = std::make_pair(x1, x2);
                for (int j = 0; j < (int)L.size(); ++j) {
                    if (L[j] == cur) {
                        tag = false;
                        break;
                    }
                }
                if (tag) L.push_back(cur);
            }
        }
        Aut = stabilized_aut;
    }

    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    for (int i = 0; i < (int)L.size(); ++i) {
        bool tag = true;
        for (int j = 0; j < (int)ordered_pairs.size(); ++j) {
            if (L[i].second == ordered_pairs[j].second) {
                tag = false;
                if (L[i].first > ordered_pairs[j].first) ordered_pairs[j].first = L[i].first;
                break;
            }
        }
        if (tag) ordered_pairs.push_back(L[i]);
    }
}

void Configuration::GraphZero_get_automorphisms(std::vector<std::vector<int>> &Aut) const {
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

double Configuration::our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int,int>> &pairs, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    int config_max_degree = get_max_degree(size, adj_mat);

    double p_size[config_max_degree];
    double pp_size[config_max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt; 
    
    p_size[0] = v_cnt;
    for (int i = 1;i < config_max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for (int i = 1; i < config_max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;
    
    std::string cur_adj_mat;
    cur_adj_mat = std::string(size*size, '0');
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int,int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    std::vector<int> tmp(size);
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i) {
            if (tmp[restricts[i].first]> tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
        }
    } while (std::next_permutation(tmp.begin(), tmp.end()));
    
    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;
    for (int i = 0; i < restricts_size; ++i) sum[i] = sum[i] /total;
    for (int i = restricts_size - 1; i> 0; --i) sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i) invariant_size[i].clear();
    
    double val = 1;
    for (int i = size - 1; i>= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)] == '1')
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)] == '1')
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j>= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)] == '1')
                invariant_size[j].push_back(c--);

        for (int j = 0; j < (int)invariant_size[i].size(); ++j)
            if (invariant_size[i][j]> 1) 
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *=  sum[j];
        if (i) val *= p_size[1] * pp_size[ cnt_forward - 1 ];
        else val *= p_size[0];
    }
    return val;
}

double Configuration::GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int,int>> &pairs, int v_cnt, unsigned int e_cnt) {
    int config_max_degree = get_max_degree(size, adj_mat);

    double p_size[config_max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;
    
    p_size[0] = v_cnt;
    for (int i = 1;i < config_max_degree; ++i) p_size[i] = p_size[i-1] * p;

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;
    
    std::string cur_adj_mat;
    cur_adj_mat = std::string(size*size, '0');
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int,int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
    
    int tmp[size];
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i) {
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) sum[i] += 1;
            else break;
        }
    } while (std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;
    
    for (int i = 0; i < restricts_size; ++i) sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i) sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i) invariant_size[i].clear();

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
        for (int j = i - 1; j >= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for (int j = 0; j < (int)invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1) 
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *=  sum[j];
        val *= p_size[cnt_forward];
    }
    return val;
}

double Configuration::Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int,int>> &pairs, int v_cnt, unsigned int e_cnt) {

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;
    
    std::string cur_adj_mat;
    cur_adj_mat = std::string(size*size, '0');
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];
    
    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());
    
    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
    int tmp[size];
    
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i) {
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) sum[i] += 1;
            else break;
        }
    } while (std::next_permutation(tmp, tmp + size));
    
    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;
    
    for (int i = 0; i < restricts_size; ++i) sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i) sum[i] /= sum[i - 1];

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *=  sum[j];
        val *= v_cnt;
        for (int j = 0; j < i - cnt_forward; ++j) val *= (1-p);
        for (int j = 0; j < cnt_forward; ++j) val *= p;
    }
    return val;
}

void Configuration::set_in_exclusion_optimize_redundancy(Graph& complete) {
    int tmp = get_in_exclusion_optimize_num();
    if (tmp <= 1) in_exclusion_optimize_redundancy = 1;
    else {
        in_exclusion_optimize_redundancy = 1;
        long long ans = pattern_matching(complete, *this);
        set_in_exclusion_optimize_num(0);
        long long true_ans = pattern_matching(complete, *this);
        set_in_exclusion_optimize_num(tmp);
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}
