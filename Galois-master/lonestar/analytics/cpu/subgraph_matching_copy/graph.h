// #pragma once
// #include "schedule.h"
// #include "vertex_set.h"
// #include "common.h"
// #include <assert.h>

// #include <cstdio>
// #include <sys/time.h>
// #include <unistd.h>
// #include <cstdlib>
// #include <omp.h>
// #include <algorithm>
// #include <cstring>
// #include <atomic>
// #include <queue>
// #include <iostream>

// #include "galois/Galois.h"
// #include "galois/AtomicHelpers.h"
// #include "galois/Reduction.h"
// #include "galois/PriorityQueue.h"
// #include "galois/Timer.h"
// #include "galois/graphs/LCGraph.h"
// #include "galois/graphs/TypeTraits.h"
// #include "Lonestar/BoilerPlate.h"
// #include "Lonestar/Utils.h"

// class Graph {
// public:
//     int v_cnt; // number of vertex
//     unsigned int e_cnt; // number of edge
//     long long tri_cnt; // number of triangle
//     double max_running_time = 60 * 60 * 24; // second

//     int *edge; // edges
//     unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]
    
//     Graph() {
//         v_cnt = 0;
//         e_cnt = 0;
//         edge = nullptr;
//         vertex = nullptr;
//     }

//     ~Graph() {
//         if(edge != nullptr) delete[] edge;
//         if(vertex != nullptr) delete[] vertex;
//     }

//     int intersection_size(int v1,int v2);
//     int intersection_size_clique(int v1,int v2);

//     //single thread triangle counting
//     long long triangle_counting();
    
//     //multi thread triangle counting
//     long long triangle_counting_mt(int thread_count);

//     //general pattern matching algorithm with multi thread
//     long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

//     //this function will be defined at code generation
//     long long unfold_pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

//     int max_degree;
// private:
//     void tc_mt(long long * global_ans);

//     void get_edge_index(int v, unsigned int& l, unsigned int& r) const;

//     void pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique = false);

//     void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

//     //this function will be defined at code generation
//     void unfold_pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

// };

// constexpr static const unsigned CHUNK_SIZE = 128U;

// int Graph::intersection_size(int v1,int v2) {
//     unsigned int l1, r1;
//     get_edge_index(v1, l1, r1);
//     unsigned int l2, r2;
//     get_edge_index(v2, l2, r2);
//     int ans = 0;
//     while(l1 < r1 && l2 < r2) {
//         if(edge[l1] < edge[l2]) {
//             ++l1;
//         }
//         else {
//             if(edge[l2] < edge[l1]) {
//                 ++l2;
//             }
//             else {
//                 // v1 and v2 have the same neighbor
//                 ++l1;
//                 ++l2;
//                 ++ans;
//             }
//         }
//     }
//     return ans;
// }

// int Graph::intersection_size_clique(int v1,int v2) {
//     unsigned int l1, r1;
//     get_edge_index(v1, l1, r1);
//     unsigned int l2, r2;
//     get_edge_index(v2, l2, r2);
//     int min_vertex = v2;
//     int ans = 0;
//     if (edge[l1] >= min_vertex || edge[l2] >= min_vertex)
//         return 0;
//     while(l1 < r1 && l2 < r2) {
//         if(edge[l1] < edge[l2]) {
//             if (edge[++l1] >= min_vertex)
//                 break;
//         }
//         else {
//             if(edge[l2] < edge[l1]) {
//                 if (edge[++l2] >= min_vertex)
//                     break;
//             }
//             else {
//                 ++ans;
//                 if (edge[++l1] >= min_vertex)
//                     break;
//                 if (edge[++l2] >= min_vertex)
//                     break;
//             }
//         }
//     }
//     return ans;
// }

// void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
// {
//     l = vertex[v];
//     r = vertex[v + 1];
// }

// long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique) {
//     long long global_ans = 0;
//     std::vector<int> allNodes;
//     for (int vertex = 0; vertex < v_cnt; ++vertex) allNodes.push_back(vertex);
//     VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
//     VertexSet subtraction_set;
//     VertexSet tmp_set;
//     subtraction_set.init();
//     long long local_ans = 0;
//     galois::do_all(
//         galois::iterate(allNodes),
//         [&](uint64_t vertex) {
//             unsigned int l, r;
//             get_edge_index(vertex, l, r);
//             for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
//             {
//                 vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id);
//             }
//             subtraction_set.push_back(vertex);
//             pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
//             subtraction_set.pop_back();
//         },
//         galois::loopname("matching"),
//         galois::chunk_size<CHUNK_SIZE>(),
//         galois::steal(),
//         galois::no_stats()
//     );
//     delete[] vertex_set;
//     // TODO : Computing multiplicty for a pattern
//     global_ans += local_ans;

//     return global_ans / schedule.get_in_exclusion_optimize_redundancy();
// }

// void Graph::pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth) // 3 same # or @ in comment are useful in code generation ###
// {
//     int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);// @@@
//     int loop_size = vertex_set[loop_set_prefix_id].get_size();
//     if (loop_size <= 0)
//         return;

//     int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
//     //Case: in_exclusion_optimize_num > 1
//     if( depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num() ) {
//         int in_exclusion_optimize_num = schedule.get_in_exclusion_optimize_num();// @@@
//         int loop_set_prefix_ids[ in_exclusion_optimize_num ];
//         loop_set_prefix_ids[0] = loop_set_prefix_id;
//         for(int i = 1; i < in_exclusion_optimize_num; ++i)
//             loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id( depth + i );
//         for(int optimize_rank = 0; optimize_rank < schedule.in_exclusion_optimize_group.size(); ++optimize_rank) {
//             const std::vector< std::vector<int> >& cur_graph = schedule.in_exclusion_optimize_group[optimize_rank];
//             long long val = schedule.in_exclusion_optimize_val[optimize_rank];
//             for(int cur_graph_rank = 0; cur_graph_rank < cur_graph.size(); ++ cur_graph_rank) {
//                 //                VertexSet tmp_set;
                
//                 //if size == 1 , we will not call intersection(...)
//                 //so we will not allocate memory for data
//                 //otherwise, we need to copy the data to do intersection(...)
//                 if(cur_graph[cur_graph_rank].size() == 1) {
//                     int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
//                     val = val * VertexSet::unorderd_subtraction_size(vertex_set[id], subtraction_set);
//                 }
//                 else {
//                     int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
//                     int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
//                     tmp_set.init(this->max_degree);
//                     tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

//                     for(int i = 2; i < cur_graph[cur_graph_rank].size(); ++i) {
//                         int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
//                         tmp_set.intersection_with(vertex_set[id]);
//                     }
//                     val = val * VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
//                 }
//                 if( val == 0 ) break;

//             }
//             local_ans += val;
//         }
//         return;// @@@
            
//     }
//     //Case: in_exclusion_optimize_num <= 1
//     if (depth == schedule.get_size() - 1)
//     {
//         // TODO : try more kinds of calculation. @@@
//         // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
//         if (schedule.get_total_restrict_num() > 0)
//         {
//             int min_vertex = v_cnt;
//             for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
//                 if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
//                     min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
//             const VertexSet& vset = vertex_set[loop_set_prefix_id];
//             int size_after_restrict = std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex) - vset.get_data_ptr();
//             if (size_after_restrict > 0)
//                 local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
//         }
//         else
//             local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set); 
//         return;// @@@
//     }
  
//     // TODO : min_vertex is also a loop invariant @@@
//     int min_vertex = v_cnt;
//     for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
//         if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
//             min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
//     int ii = 0;
//     for (int &i = ii; i < loop_size; ++i)
//     {
//         if (min_vertex <= loop_data_ptr[i])
//             break;
//         int vertex = loop_data_ptr[i];
//         if (subtraction_set.has_data(vertex))
//             continue;
//         unsigned int l, r;
//         get_edge_index(vertex, l, r);
//         bool is_zero = false;
//         for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
//         {
//             vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, vertex);
//             if( vertex_set[prefix_id].get_size() == 0) {
//                 is_zero = true;
//                 break;
//             }
//         }
//         if( is_zero ) continue;
//         //subtraction_set.insert_ans_sort(vertex);
//         subtraction_set.push_back(vertex);
//         pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
//         subtraction_set.pop_back(); // @@@
//     }
// } 
// // ###
